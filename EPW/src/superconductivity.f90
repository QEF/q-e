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
    !-----------------------------------------------------------------------
    SUBROUTINE eliashberg_iso_iaxis
    !-----------------------------------------------------------------------
    !!
    !! This routine is the driver of the self-consistent cycle for the isotropic 
    !! Eliashberg equations on the imaginary-axis.  
    !!
    !
    USE kinds,         ONLY : DP
    USE io_global,     ONLY : stdout
    USE control_flags, ONLY : iverbosity
    USE epwcom,        ONLY : nsiter, nstemp, broyden_beta, broyden_ndim, &
                              limag, lpade, lacon
    USE eliashbergcom, ONLY : nsw, nsiw, Deltai, Deltaip, Delta, Deltap, estemp
    USE constants_epw, ONLY : kelvin2eV, ci
    USE mp,            ONLY : mp_bcast, mp_barrier, mp_sum
    ! 
    IMPLICIT NONE
    !
    ! Local variables
    INTEGER :: itemp
    !! Counter on temperature index
    INTEGER :: iter
    !! Counter on iteration steps
    INTEGER :: N
    !! Maximum frequency 
    REAL(DP) :: tcpu, rdeltaout(nsw), rdeltain(nsw), cdeltaout(nsw), cdeltain(nsw)
    REAL(DP), EXTERNAL :: get_clock
    LOGICAL :: conv
    !
    CALL start_clock( 'iso_iaxis' )
    !
    DO itemp = 1, nstemp ! loop over temperature
      !
      WRITE(stdout,'(a)') '    '
      WRITE(stdout,'(5x,a,i3,a,f12.5,a,a,i3,a)') 'temp(', itemp, ') = ', estemp(itemp)/kelvin2eV, ' K '
      WRITE(stdout,'(a)') '    '
      WRITE(stdout,'(5x,a)') 'Solve isotropic Eliashberg equations on imaginary-axis' 
      WRITE(stdout,'(a)') '    '
      WRITE(stdout,'(5x,a,i6,a,i6)') 'Total number of frequency points nsiw ( ', itemp, ' ) = ', nsiw(itemp)
      WRITE(stdout,'(a)') '    '
      CALL start_clock( 'iaxis_imag' )
      CALL gen_freqgrid_iaxis( itemp )
      !
      IF ( limag ) THEN
        iter = 1
        conv = .false.
        DO WHILE ( .not. conv .AND. iter .le. nsiter )
          CALL sum_eliashberg_iso_iaxis( itemp, iter, conv )
          CALL mix_broyden( nsiw(itemp), Deltai, Deltaip, broyden_beta, iter, broyden_ndim, conv )
          iter = iter + 1
        ENDDO ! iter
        !
        IF ( conv ) THEN
          !
          ! SP : Only print the Free energy if the user want it
          !
          IF ( iverbosity .eq. 2 ) THEN
            CALL free_energy( itemp )
          ENDIF
          WRITE(stdout,'(a)') '  '
          CALL stop_clock( 'iaxis_imag' )
          CALL print_clock( 'iaxis_imag' )
        ELSEIF ( .not. conv .AND. (iter-1) .eq. nsiter ) THEN
          CALL deallocate_eliashberg
          WRITE(stdout,'(a)') 'not converged  '
          CALL stop_clock( 'iaxis_imag' )
          CALL print_clock( 'iaxis_imag' )
          CALL errore('sum_eliashberg_iso_iaxis','converged was not reached',1)
          RETURN
        ENDIF
      ENDIF
      !
      IF ( lpade ) THEN 
        WRITE(stdout,'(a)') '    '
        WRITE(stdout,'(5x,a)') 'Pade approximant of isotropic Eliashberg equations from imaginary-axis to real-axis'
        WRITE(stdout,'(a)') '    '
        CALL start_clock( 'raxis_pade' )
        !
        iter = 1
        conv = .false.
        N = 80 * nsiw(itemp) / 100
        IF ( mod(N,2) .ne. 0 ) N = N + 1
        DO WHILE ( .not. conv .AND. iter .le. nsiter )
          CALL pade_cont_iso_iaxis_to_raxis( itemp, N, conv )
          N = N - 2
          iter = iter + 1
        ENDDO
        !
        IF ( conv ) THEN
          CALL dos_quasiparticle( itemp )
          WRITE(stdout,'(a)') '  '
          CALL stop_clock( 'raxis_pade' )
          CALL print_clock( 'raxis_pade' )
          WRITE(stdout,'(a)') '  '
        ELSEIF ( .not. conv  .AND. (iter-1) .eq. nsiter ) THEN
          CALL deallocate_eliashberg
          WRITE(stdout,'(a)') '  '
          CALL stop_clock( 'raxis_pade' )
          CALL print_clock( 'raxis_pade' )
          CALL errore('pade_cont_iso_iaxis_to_raxis','converged was not reached',1)
          RETURN
        ENDIF
      ENDIF 
      !
      IF ( lacon ) THEN 
        WRITE(stdout,'(a)') '    '
        WRITE(stdout,'(5x,a)') 'Analytic continuation of isotropic Eliashberg equations from imaginary-axis to real-axis'
        WRITE(stdout,'(a)') '    '
        WRITE(stdout,'(5x,a,i6)') 'Total number of frequency points nsw = ', nsw
        WRITE(stdout,'(a)') '    '
        CALL start_clock( 'raxis_acon' )
        !
        iter = 1
        conv = .false.
        DO WHILE ( .not. conv .AND. iter .le. nsiter )
          CALL analytic_cont_iso_iaxis_to_raxis( itemp, iter, conv )
          rdeltain(:)  = real(Deltap(:))
          cdeltain(:)  = aimag(Deltap(:))
          rdeltaout(:) = real(Delta(:))
          cdeltaout(:) = aimag(Delta(:))
          CALL mix_broyden ( nsw, rdeltaout, rdeltain, broyden_beta, iter, broyden_ndim, conv )
          CALL mix_broyden2( nsw, cdeltaout, cdeltain, broyden_beta, iter, broyden_ndim, conv )
          Deltap(:) = rdeltain(:) + ci * cdeltain(:)
          iter = iter + 1
        ENDDO ! iter
        !
        IF ( conv ) THEN
          CALL dos_quasiparticle( itemp )
          WRITE(stdout,'(a)') ' '
          CALL stop_clock( 'raxis_acon' )
          CALL print_clock( 'raxis_acon' )
          WRITE(stdout,'(a)') ' '
        ELSEIF ( .not. conv .AND. (iter-1) .eq. nsiter ) THEN
          CALL deallocate_eliashberg
          WRITE(stdout,'(a)') '  '
          CALL stop_clock( 'raxis_acon' )
          CALL print_clock( 'raxis_acon' )
          CALL errore('analytic_cont_iso_iaxis_to_raxis','converged was not reached',1)
          RETURN
        ENDIF
      ENDIF
      !
      CALL deallocate_eliashberg_iso_iaxis
      !
      tcpu = get_clock('iso_iaxis')
      WRITE(stdout,'(5x,a,i3,a,f8.1,a)') 'itemp = ', itemp, '   total cpu time :', tcpu, ' secs'
      !
    ENDDO ! itemp
    !
    CALL stop_clock( 'iso_iaxis' )
    !
    RETURN
    !
    END SUBROUTINE eliashberg_iso_iaxis
    !
    !-----------------------------------------------------------------------
    SUBROUTINE sum_eliashberg_iso_iaxis( itemp, iter, conv ) 
    !-----------------------------------------------------------------------
    !!
    !! This routine solves the isotropic Eliashberg equations on the imaginary-axis
    !!
    !
    USE kinds,         ONLY : DP
    USE io_global,     ONLY : stdout
    USE epwcom,        ONLY : nsiter, nstemp, muc, conv_thr_iaxis
    USE eliashbergcom, ONLY : nsiw, estemp, gap0, gap, wsi, NZnormi, Znormi, Deltai, Deltaip, Keri
    USE constants_epw, ONLY : pi
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
    REAL(DP) :: esqrt, kernelp, kernelm, lambdap, lambdam, absdelta, reldelta, errdelta
    REAL(DP), ALLOCATABLE :: wesqrt(:), desqrt(:)
    REAL(DP), ALLOCATABLE, SAVE :: Deltaold(:)
    !
    IF ( .not. ALLOCATED(wesqrt) ) ALLOCATE( wesqrt(nsiw(itemp)) )
    IF ( .not. ALLOCATED(desqrt) ) ALLOCATE( desqrt(nsiw(itemp)) )
    !
    IF ( iter .eq. 1 ) THEN
      IF ( .not. ALLOCATED(gap) )      ALLOCATE( gap(nstemp) )
      IF ( .not. ALLOCATED(Deltai) )   ALLOCATE( Deltai(nsiw(itemp)) )
      IF ( .not. ALLOCATED(Deltaip) )  ALLOCATE( Deltaip(nsiw(itemp)) )
      IF ( .not. ALLOCATED(Znormi) )   ALLOCATE( Znormi(nsiw(itemp)) )
      IF ( .not. ALLOCATED(NZnormi) )  ALLOCATE( NZnormi(nsiw(itemp)) )
      gap(itemp) = 0.d0
      Deltaip(:) = 0.d0
      Deltaip(:) = gap0
      !
      CALL kernel_iso_iaxis( itemp )
    ENDIF
    Znormi(:) = 0.d0
    NZnormi(:) = 0.d0
    Deltai(:) = 0.d0
    !
    IF ( iter .eq. 1 ) THEN
      IF ( .not. ALLOCATED(Deltaold) ) ALLOCATE( Deltaold(nsiw(itemp)) )
      Deltaold(:) = gap0
    ENDIF
    absdelta = 0.d0   
    reldelta = 0.d0 
    DO iw = 1, nsiw(itemp) ! loop over omega
      DO iwp = 1, nsiw(itemp) ! loop over omega_prime
        ! this step is performed at each iter step only for iw=1 since it is independ of wsi(iw)
        IF ( iw .eq. 1 ) THEN
          esqrt = 1.d0 / sqrt( wsi(iwp)**2.d0 + Deltaip(iwp)**2.d0 )
          wesqrt(iwp) =  wsi(iwp) * esqrt 
          desqrt(iwp) =  Deltaip(iwp) * esqrt 
        ENDIF
        lambdam = Keri( abs(iw-iwp)+1 )
        lambdap = Keri( abs(iw+iwp) )
        kernelm = lambdam - lambdap
        kernelp = lambdam + lambdap
        NZnormi(iw) = NZnormi(iw) + kernelm
        Znormi(iw) = Znormi(iw) + wesqrt(iwp) * kernelm 
        Deltai(iw) = Deltai(iw) + desqrt(iwp) * ( kernelp - 2.d0 * muc ) 
      ENDDO ! iwp
      Znormi(iw) = 1.d0 + pi * estemp(itemp) * Znormi(iw) / wsi(iw)
      NZnormi(iw) = 1.d0 + pi * estemp(itemp) * NZnormi(iw) / wsi(iw)
      Deltai(iw) = pi * estemp(itemp) * Deltai(iw) / Znormi(iw)
      reldelta = reldelta + abs( Deltai(iw) - Deltaold(iw) )
      absdelta = absdelta + abs( Deltai(iw) ) 
    ENDDO ! iw 
    errdelta = reldelta / absdelta
    Deltaold(:) = Deltai(:)
    !
    WRITE(stdout,'(5x,a,i6,a,ES20.10,a,ES20.10,a,ES20.10)') 'iter = ', iter, '   error = ', errdelta, &
                                           '   Znormi(1) = ', Znormi(1), '   Deltai(1) = ', Deltai(1)
    !
    IF ( errdelta .lt. conv_thr_iaxis ) conv = .true.
    IF ( errdelta .lt. conv_thr_iaxis .OR. iter .eq. nsiter ) THEN
      gap(itemp) = Deltai(1)
      gap0 = gap(itemp)
      CALL eliashberg_write_iaxis( itemp )
    ENDIF
    !
    IF( ALLOCATED(wesqrt) ) DEALLOCATE(wesqrt)
    IF( ALLOCATED(desqrt) ) DEALLOCATE(desqrt)
    !
    IF ( conv .OR. iter .eq. nsiter ) THEN
      IF( ALLOCATED(Deltaold) ) DEALLOCATE(Deltaold)
      WRITE(stdout,'(5x,a,i6)') 'Convergence was reached in nsiter = ', iter
    ENDIF
    IF ( .not. conv .AND. iter .eq. nsiter ) THEN
      WRITE(stdout,'(5x,a,i6)') 'Convergence was not reached in nsiter = ', iter
      CALL errore('sum_eliashberg_iso_iaxis','increase nsiter or reduce conv_thr_iaxis',1)
      !CALL deallocate_eliashberg
    ENDIF
    !
    RETURN
    !
    END SUBROUTINE sum_eliashberg_iso_iaxis
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
    !
    !-----------------------------------------------------------------------
    SUBROUTINE analytic_cont_iso_iaxis_to_raxis( itemp, iter, conv ) 
    !-----------------------------------------------------------------------
    !!
    !! This routine does the analyic continuation of the isotropic Eliashberg equations 
    !! from the imaginary-axis to the real axis
    !! reference F. Marsiglio, M. Schossmann, and J. Carbotte, Phys. Rev. B 37, 4965 (1988)
    !!
    !
    USE kinds,         ONLY : DP
    USE io_global,     ONLY : stdout
    USE epwcom,        ONLY : nqstep, nsiter, conv_thr_racon, lpade
    USE eliashbergcom, ONLY : nsw, estemp, dwsph, ws, gap, a2f_iso, Dsumi, Zsumi, & 
                              Delta, Deltap, Znorm, Znormp, Gp, Gm
    USE constants_epw, ONLY : pi, ci
    ! 
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: itemp
    !! Counter on iteration
    INTEGER, INTENT(in) :: iter
    !! Counter on the iteration number
    LOGICAL, INTENT(inout) :: conv
    !! True if the calculation is converged
    ! 
    ! Local variables
    INTEGER :: i, iw, iwp 
    REAL(kind=DP) :: rgammap, rgammam, absdelta, reldelta, errdelta
    COMPLEX(kind=DP) :: esqrt, root
    COMPLEX(DP), ALLOCATABLE, SAVE :: Deltaold(:)
    CHARACTER (len=256) :: cname
    !
    IF ( iter .eq. 1 ) THEN
      IF ( .not. ALLOCATED(Delta) )    ALLOCATE( Delta(nsw) )
      IF ( .not. ALLOCATED(Deltap) )   ALLOCATE( Deltap(nsw) )
      IF ( .not. ALLOCATED(Znorm) )    ALLOCATE( Znorm(nsw) )
      IF ( .not. ALLOCATED(Znormp) )   ALLOCATE( Znormp(nsw) )
      IF ( .not. ALLOCATED(Deltaold) ) ALLOCATE( Deltaold(nsw) )
      Deltap(:) = (0.d0, 0.d0)
      Deltaold(:) = (0.d0, 0.d0)
      IF ( lpade ) THEN
         Deltap(:) = Delta(:)
         Deltaold(:) = Delta(:)
      ELSE 
         Deltap(:) = gap(itemp)
         Deltaold(:) = gap(itemp)
      ENDIF
      Znormp(:) = (1.d0, 0.d0)
      IF ( .not. ALLOCATED(Gp) ) ALLOCATE( Gp(nsw,nqstep) )
      IF ( .not. ALLOCATED(Gm) ) ALLOCATE( Gm(nsw,nqstep) )
      IF ( .not. ALLOCATED(Dsumi) ) ALLOCATE( Dsumi(nsw) )
      IF ( .not. ALLOCATED(Zsumi) ) ALLOCATE( Zsumi(nsw) )
      !CALL kernel_iso_iaxis_analytic_cont( itemp, Zsumi, Dsumi )
      !SP -debug
      CALL kernel_iso_iaxis_analytic_cont( itemp )
    ENDIF
    Znorm(:) = (0.d0, 0.d0)
    Delta(:) = (0.d0, 0.d0)
    !
    absdelta = 0.d0   
    reldelta = 0.d0
    DO iw = 1, nsw ! loop over omega
      DO iwp = 1, nqstep ! loop over omega_prime
        IF ( iter .eq. 1 ) THEN 
          CALL gamma_acont( ws(iw), ws(iwp), estemp(itemp), rgammap, rgammam )
          Gp(iw,iwp) = rgammap
          Gm(iw,iwp) = rgammam
        ENDIF
        !
        i = iw + iwp - 1
        IF ( i .le. nsw ) THEN
          root = sqrt( Znormp(i)**2.d0 * ( ws(i)**2.d0 - Deltap(i)**2.d0 ) )
          IF ( aimag(root) .lt. 0.d0 ) THEN 
             esqrt = Znormp(i) / conjg(root)
          ELSE  
             esqrt = Znormp(i) / root
          ENDIF
          esqrt = esqrt * Gp(iw,iwp) * a2f_iso(iwp) 
          Znorm(iw) = Znorm(iw) - ws(i) * esqrt 
          Delta(iw) = Delta(iw) - Deltap(i) * esqrt 
        ENDIF
        ! 
        i = abs(iw - iwp) + 1
        root = sqrt( Znormp(i)**2.d0 * ( ws(i)**2.d0 - Deltap(i)**2.d0 ) )
        IF ( aimag(root) .lt. 0.d0 ) THEN 
           esqrt = Znormp(i) / conjg(root)
        ELSE  
           esqrt = Znormp(i) / root
        ENDIF
        esqrt = esqrt * Gm(iw,iwp) * a2f_iso(iwp) 
        IF ( iw .lt. iwp ) THEN 
           Znorm(iw) = Znorm(iw) - ws(i) * esqrt 
        ELSE
           Znorm(iw) = Znorm(iw) + ws(i) * esqrt 
        ENDIF
        Delta(iw) = Delta(iw) + Deltap(i) * esqrt
      ENDDO ! iwp
      Znorm(iw) = 1.d0 + pi * ( - estemp(itemp) * Zsumi(iw) + ci * Znorm(iw) * dwsph ) / ws(iw)
      Delta(iw) = pi * ( estemp(itemp) * Dsumi(iw) + ci * Delta(iw) * dwsph ) / Znorm(iw)
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
       IF( ALLOCATED(Deltaold) ) DEALLOCATE(Deltaold)
       WRITE(stdout,'(5x,a,i6)') 'Convergence was reached in nsiter = ', iter
    ENDIF
    IF ( .not. conv .AND. iter .eq. nsiter ) THEN
       WRITE(stdout,'(5x,a,i6)') 'Convergence was not reached in nsiter = ', iter
       CALL errore('analytic_cont_iso_iaxis_to_raxis','increase nsiter or reduce conv_thr_racon',-1)
    ENDIF
    !
    RETURN
    !
    END SUBROUTINE analytic_cont_iso_iaxis_to_raxis
    !
    !-----------------------------------------------------------------------
    SUBROUTINE pade_cont_iso_iaxis_to_raxis( itemp, N, conv )
    !-----------------------------------------------------------------------
    !
    ! This routine uses pade approximants to continue the isotropic Eliashberg equations 
    ! from the imaginary-axis to the real-axis
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
    USE eliashbergcom, ONLY : nsw, ws, wsi, gap, Delta, Znorm, Deltai, Znormi
    USE constants_epw, ONLY : cone, ci
    ! 
    IMPLICIT NONE
    !
    INTEGER, INTENT (in) :: itemp
    !! Counter on temperature
    INTEGER, INTENT (in) :: N
    !! Maximum number of frequency 
    LOGICAL, INTENT (inout) :: conv
    !! True if the calculation is converged
    ! 
    ! Local variable
    INTEGER :: iw
    REAL(DP) :: absdelta, reldelta, errdelta
    COMPLEX(DP) :: a(N), b(N), z(N), u(N), v(N)
    COMPLEX(DP) :: omega, padapp, Deltaold(nsw)
    CHARACTER (len=256) :: cname
    !
    Deltaold(:) = gap(itemp)
    absdelta = 0.d0
    reldelta = 0.d0
    !
    IF ( .not. ALLOCATED(Delta) )  ALLOCATE( Delta(nsw) )
    IF ( .not. ALLOCATED(Znorm) )  ALLOCATE( Znorm(nsw) )
    Znorm(:) = (0.d0, 0.d0)
    Delta(:) = (0.d0, 0.d0)
    !
    DO iw = 1, N
       z(iw) = ci * wsi(iw)
       u(iw) = cone * Deltai(iw)
       v(iw) = cone * Znormi(iw)
    ENDDO
    !
    CALL pade_coeff( N, z, u, a )
    CALL pade_coeff( N, z, v, b )
    !
    DO iw = 1, nsw
       omega = cone * ws(iw)
       CALL pade_eval( N, z, a, omega, padapp )
       Delta(iw) = padapp
       CALL pade_eval( N, z, b, omega, padapp )
       Znorm(iw) = padapp
       reldelta = reldelta + abs( Delta(iw) - Deltaold(iw) ) 
       absdelta = absdelta + abs( Delta(iw) )
    ENDDO
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
    !
  !  IF ( .not. conv ) THEN
  !     WRITE(stdout,'(5x,a,i6)') 'Convergence was not reached pade = ', N
  !     CALL errore('pade_cont_iso_iaxis_to_raxis','decrease number of Pade approximants',-1)
  !  ENDIF
    !
    RETURN
    !
    END SUBROUTINE pade_cont_iso_iaxis_to_raxis
    !
    !-----------------------------------------------------------------------
    !                                                                            
    !-----------------------------------------------------------------------
    SUBROUTINE eliashberg_aniso_iaxis
    !-----------------------------------------------------------------------
    !
    ! This routine is the driver of the self-consistent cycle for the anisotropic 
    ! Eliashberg equations on the imaginary-axis.  
    !
    USE kinds,         ONLY : DP
    USE io_global,     ONLY : stdout
    USE control_flags, ONLY : iverbosity
    USE epwcom,        ONLY : nsiter, nstemp, broyden_beta, broyden_ndim, & 
                              limag, lpade, lacon, fsthick, imag_read, wscut
    USE eliashbergcom, ONLY : nsw, nsiw, ADelta, ADeltap, ADeltai, ADeltaip, &
                              estemp, nkfs, nbndfs, ekfs, ef0
    USE constants_epw, ONLY : kelvin2eV, ci, pi
    USE io_global,     ONLY : ionode_id
    USE mp_global,     ONLY : inter_pool_comm
    USE mp,            ONLY : mp_bcast, mp_barrier
    USE mp_world,      ONLY : mpime
    ! 
    IMPLICIT NONE
    !
    INTEGER :: itemp, iter, N, ik, ibnd, imelt
    REAL(DP) :: tcpu, rdeltaout(nsw), rdeltain(nsw), cdeltaout(nsw), cdeltain(nsw)
    REAL(DP), EXTERNAL :: get_clock
    LOGICAL :: conv 
    !
    CALL start_clock( 'aniso_iaxis' )
    !
    DO itemp = 1, nstemp ! loop over temperature
       !
       WRITE(stdout,'(a)') '  '
       WRITE(stdout,'(5x,a,i3,a,f8.4,a,a,i3,a)') 'temp(', itemp, ') = ', estemp(itemp)/kelvin2eV, ' K '
       WRITE(stdout,'(a)') '  '
       IF ( limag .AND. .not. imag_read ) THEN 
          WRITE(stdout,'(5x,a)') 'Solve anisotropic Eliashberg equations on imaginary-axis ' 
       ELSEIF ( limag .AND. imag_read ) THEN
          WRITE(stdout,'(5x,a)') 'Read from file Delta and Znorm on imaginary-axis '
       ENDIF
       WRITE(stdout,'(a)') '  '
       WRITE(stdout,'(5x,a,i6,a,i6)') 'Total number of frequency points nsiw ( ', itemp, ' ) = ', nsiw(itemp)
       WRITE(stdout,'(5x,a,f10.4)') 'Cutoff frequency wscut = ', (2.d0*nsiw(itemp)+1)*pi*estemp(itemp)
       WRITE(stdout,'(a)') '  '
       CALL start_clock( 'iaxis_imag' )
       CALL gen_freqgrid_iaxis( itemp )
       !
       IF ( ( limag .AND. .not. imag_read ) .OR. ( limag .AND. imag_read .AND. itemp .ne. 1 ) ) THEN
          iter = 1
          conv = .false.
          DO WHILE ( .not. conv .AND. iter .le. nsiter ) 
             CALL sum_eliashberg_aniso_iaxis( itemp, iter, conv )
             IF (mpime .eq. ionode_id) THEN
                DO ik = 1, nkfs
                   DO ibnd = 1, nbndfs
                      IF ( abs( ekfs(ibnd,ik) - ef0 ) .lt. fsthick ) THEN
                         CALL mix_broyden_aniso( ik, ibnd, nsiw(itemp), & 
                              ADeltai(ibnd,ik,:), ADeltaip(ibnd,ik,:), broyden_beta, iter, broyden_ndim, conv )
                      ENDIF
                   ENDDO
                ENDDO
             ENDIF
             CALL mp_bcast( ADeltai, ionode_id, inter_pool_comm )
             CALL mp_bcast( ADeltaip, ionode_id, inter_pool_comm )
             CALL mp_barrier(inter_pool_comm)
             iter = iter + 1
          ENDDO ! iter
          !
          IF ( conv ) THEN
             IF ( ALLOCATED(ADeltaip) ) DEALLOCATE(ADeltaip)
             !
             ! SP : Only print the Free energy if the user want it
             !
             IF ( iverbosity .eq. 2 ) THEN
                IF (mpime .eq. ionode_id) THEN
                   CALL free_energy( itemp )
                ENDIF
                CALL mp_barrier(inter_pool_comm)
             ENDIF
             !
             WRITE(stdout,'(a)') '  '
             CALL stop_clock( 'iaxis_imag' )
             CALL print_clock( 'iaxis_imag' )
          ELSEIF ( .not. conv .AND. (iter-1) .eq. nsiter ) THEN
             CALL deallocate_eliashberg
             WRITE(stdout,'(a)') 'not converged  '
             CALL stop_clock( 'iaxis_imag' )
             CALL print_clock( 'iaxis_imag' )
             CALL errore('sum_eliashberg_aniso_iaxis','convergence was not reached',1)
             RETURN
          ENDIF
       ELSEIF ( limag .AND. imag_read .AND. itemp .eq. 1 ) THEN
          CALL eliashberg_read_aniso_iaxis( itemp )
       ENDIF
       !
       IF ( lpade ) THEN 
          WRITE(stdout,'(a)') '  '  
          WRITE(stdout,'(5x,a)') 'Pade approximant of anisotropic Eliashberg equations from imaginary-axis to real-axis'
          WRITE(stdout,'(5x,a,f10.4)') 'Cutoff frequency wscut = ', wscut
          WRITE(stdout,'(a)') '  '
          CALL start_clock( 'raxis_pade' )
          conv = .false.
          N = 90 * nsiw(itemp) / 100
          IF ( mod(N,2) .ne. 0 ) N = N + 1
          CALL pade_cont_aniso_iaxis_to_raxis( itemp, N, conv )
          !
          IF ( conv ) THEN
             IF (mpime .eq. ionode_id) THEN
                CALL dos_quasiparticle( itemp )
             ENDIF
             CALL mp_barrier(inter_pool_comm)
             CALL stop_clock( 'raxis_pade' )
             CALL print_clock( 'raxis_pade' )
             WRITE(stdout,'(a)') '  '
          ELSEIF ( .not. conv  .AND. (iter-1) .eq. nsiter ) THEN
             CALL deallocate_eliashberg
             WRITE(stdout,'(a)') '  '
             CALL stop_clock( 'raxis_pade' )
             CALL print_clock( 'raxis_pade' )
             CALL errore('pade_cont_iso_iaxis_to_raxis','converged was not reached',1)
             RETURN
          ENDIF
       ENDIF
       !
       IF ( lacon ) THEN 
          WRITE(stdout,'(a)') '  '
          WRITE(stdout,'(5x,a)') 'Analytic continuation of anisotropic Eliashberg equations from imaginary-axis to real-axis'
          WRITE(stdout,'(a)') '  '
          WRITE(stdout,'(5x,a,i6)') 'Total number of frequency points nsw = ', nsw
          WRITE(stdout,'(5x,a,f10.4)') 'Cutoff frequency wscut = ', wscut
          WRITE(stdout,'(a)') '    '
          CALL start_clock( 'raxis_acon' )
          !
          iter = 1
          conv = .false.
          DO WHILE ( .not. conv .AND. iter .le. nsiter )
            CALL analytic_cont_aniso_iaxis_to_raxis( itemp, iter, conv )
            IF (mpime .eq. ionode_id) THEN
              DO ik = 1, nkfs
                DO ibnd = 1, nbndfs
                  IF ( abs( ekfs(ibnd,ik) - ef0 ) .lt. fsthick ) THEN
                    rdeltain(:)  = real(ADeltap(ibnd,ik,:))
                    cdeltain(:)  = aimag(ADeltap(ibnd,ik,:))
                    rdeltaout(:) = real(ADelta(ibnd,ik,:))
                    cdeltaout(:) = aimag(ADelta(ibnd,ik,:))
                    CALL mix_broyden_aniso ( ik, ibnd, nsw, rdeltaout, rdeltain, broyden_beta, iter, broyden_ndim, conv )
                    CALL mix_broyden2_aniso( ik, ibnd, nsw, cdeltaout, cdeltain, broyden_beta, iter, broyden_ndim, conv )
                    ADeltap(ik,ibnd,:) = rdeltain(:) + ci * cdeltain(:)
                  ENDIF
                ENDDO
              ENDDO
            ENDIF
            CALL mp_bcast( ADelta, ionode_id, inter_pool_comm )
            CALL mp_bcast( ADeltap, ionode_id, inter_pool_comm )
            CALL mp_barrier(inter_pool_comm)
            iter = iter + 1
          ENDDO ! iter
          !
          IF ( conv ) THEN 
            IF (mpime .eq. ionode_id) THEN
               CALL dos_quasiparticle( itemp )
            ENDIF
            CALL mp_barrier(inter_pool_comm)
            WRITE(stdout,'(a)') '  '
            CALL stop_clock( 'raxis_acon' )
            CALL print_clock( 'raxis_acon' )
          ELSEIF ( .not. conv .AND. (iter-1) .eq. nsiter ) THEN
            CALL deallocate_eliashberg
            WRITE(stdout,'(a)') '  '
            CALL stop_clock( 'raxis_acon' )
            CALL print_clock( 'raxis_acon' )
            CALL errore('analytic_cont_aniso_iaxis_to_raxis','convergence was not reached',1)
            RETURN
          ENDIF
          !
       ENDIF
       !
       CALL deallocate_eliashberg_aniso_iaxis
       ! remove memory allocated for wsi, Deltai, Znormi, NZnormi, ADeltai, AZnormi, NAZnormi
       imelt = ( 4 + 3 * nbndfs * nkfs ) * nsiw(itemp) 
       CALL mem_size_eliashberg( -imelt )
       !
       CALL deallocate_eliashberg_aniso_raxis
       IF ( lpade ) THEN
         ! remove memory allocated for ws, Delta, Znorm, ADelta, AZnorm
         imelt = nsw + 2 * ( 2 + 2 * nbndfs * nkfs ) * nsw
         CALL mem_size_eliashberg( -imelt )
       ELSEIF ( lacon ) THEN 
         ! remove memory allocated for ws, Delta, Znorm, ADelta, ADeltap, AZnorm, AZnormp
         imelt = nsw + 2 * ( 2 + 4 * nbndfs * nkfs ) * nsw
         CALL mem_size_eliashberg( -imelt )
       ENDIF
       ! 
       tcpu = get_clock('aniso_iaxis')
       WRITE(stdout,'(5x,a,i3,a,f18.2,a)') 'itemp = ', itemp, '   total cpu time :', tcpu, ' secs'
       !
    ENDDO ! itemp
    !
    CALL stop_clock( 'aniso_iaxis' )
    !
    RETURN
    !
    END SUBROUTINE eliashberg_aniso_iaxis
    !
    !-----------------------------------------------------------------------
    SUBROUTINE sum_eliashberg_aniso_iaxis( itemp, iter, conv ) 
    !-----------------------------------------------------------------------
    !
    ! This routine solves the anisotropic Eliashberg equations on the imaginary-axis
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
    USE elph2,         ONLY : wqf
    USE epwcom,        ONLY : nsiter, nstemp, muc, conv_thr_iaxis, fsthick
    USE eliashbergcom, ONLY : nsiw, estemp, gap0, gap, Agap, wsi, AKeri, limag_fly, & 
                              NAZnormi, AZnormi, ADeltai, ADeltaip, NZnormi, Znormi, & 
                              Deltai, wsphmax, nkfs, nbndfs, dosef, ef0, ixkqf, ixqfs, & 
                              nqfs, wkfs, w0g, ekfs
    USE constants_epw, ONLY : pi  
    USE io_global,     ONLY : ionode_id
    USE mp_global,     ONLY : inter_pool_comm
    USE mp_world,      ONLY : mpime
    USE mp,            ONLY : mp_bcast, mp_barrier, mp_sum
    ! 
    IMPLICIT NONE
    !
    INTEGER  :: iw, iwp, itemp, iter, ik, iq, iq0, ibnd, jbnd, & 
               lower_bnd, upper_bnd, imelt
    REAL(DP) :: esqrt, absdelta, reldelta, errdelta, weight
    REAL(DP) :: kernelp, kernelm, lambdap, lambdam
    REAL(DP), ALLOCATABLE :: wesqrt(:,:,:), desqrt(:,:,:)
    REAL(DP), ALLOCATABLE, SAVE :: Deltaold(:)
    LOGICAL  :: conv
    !
    IF ( .not. ALLOCATED(wesqrt) ) ALLOCATE( wesqrt(nbndfs,nkfs,nsiw(itemp)) )
    IF ( .not. ALLOCATED(desqrt) ) ALLOCATE( desqrt(nbndfs,nkfs,nsiw(itemp)) )
    wesqrt(:,:,:) = 0.d0
    desqrt(:,:,:) = 0.d0
    !
    IF ( iter .eq. 1 ) THEN
       !
       IF ( itemp .eq. 1 ) THEN 
          ! get the size of required memory for  gap, Agap
          imelt = ( 1 + nbndfs * nkfs ) * nstemp 
          CALL mem_size_eliashberg( imelt )
       ENDIF
       !
       ! get the size of required memory for  
       ! wesqrt, desqrt, Deltai, Znormi, NZnormi, ADeltai, ADeltaip, AZnormi, NAZnormi, Deltaold
       imelt = ( 4 + 6 * nbndfs * nkfs ) * nsiw(itemp)
       CALL mem_size_eliashberg( imelt )
       !
       IF ( .not. ALLOCATED(gap) )       ALLOCATE( gap(nstemp) )
       IF ( .not. ALLOCATED(Agap) )      ALLOCATE( Agap(nbndfs,nkfs,nstemp) )
       IF ( .not. ALLOCATED(Deltai) )    ALLOCATE( Deltai(nsiw(itemp)) )
       IF ( .not. ALLOCATED(Znormi) )    ALLOCATE( Znormi(nsiw(itemp)) )
       IF ( .not. ALLOCATED(NZnormi) )   ALLOCATE( NZnormi(nsiw(itemp)) )
       IF ( .not. ALLOCATED(ADeltai) )   ALLOCATE( ADeltai(nbndfs,nkfs,nsiw(itemp)) )
       IF ( .not. ALLOCATED(ADeltaip) )  ALLOCATE( ADeltaip(nbndfs,nkfs,nsiw(itemp)) )
       IF ( .not. ALLOCATED(AZnormi) )   ALLOCATE( AZnormi(nbndfs,nkfs,nsiw(itemp)) )
       IF ( .not. ALLOCATED(NAZnormi) )  ALLOCATE( NAZnormi(nbndfs,nkfs,nsiw(itemp)) )
       gap(itemp) = 0.d0
       Agap(:,:,itemp) = 0.d0
       ADeltaip(:,:,:) = 0.d0
       !
       DO ik = 1, nkfs
          DO ibnd = 1, nbndfs
             IF ( abs( ekfs(ibnd,ik) - ef0 ) .lt. fsthick ) THEN
                DO iw = 1, nsiw(itemp)
                   IF ( wsi(iw) .lt. 2.d0*wsphmax ) THEN
                      ADeltaip(ibnd,ik,iw) = gap0
                   ELSE
                      ADeltaip(ibnd,ik,iw) = 0.d0
                   ENDIF
                ENDDO
             ENDIF
          ENDDO ! ibnd
       ENDDO ! ik
       !
       CALL eliashberg_memlt_aniso_iaxis( itemp )
       IF ( .not. limag_fly ) CALL kernel_aniso_iaxis( itemp )
       !
    ENDIF 
    Deltai(:) = 0.d0
    Znormi(:) = 0.d0
    NZnormi(:) = 0.d0
    ADeltai(:,:,:) = 0.d0
    AZnormi(:,:,:) = 0.d0
    NAZnormi(:,:,:) = 0.d0
    !
    CALL fkbounds( nkfs, lower_bnd, upper_bnd )
    !
    DO ik = lower_bnd, upper_bnd
       DO ibnd = 1, nbndfs
          IF ( ( abs( ekfs(ibnd,ik) - ef0 ) .lt. fsthick ) ) THEN
             DO iq = 1, nqfs(ik)
                ! iq0 - index of q-point on the full q-mesh
                iq0 = ixqfs(ik,iq)
                DO jbnd = 1, nbndfs
                   IF ( ( abs( ekfs(jbnd,ixkqf(ik,iq0)) - ef0 ) .lt. fsthick ) ) THEN
                      weight = wqf(iq) * w0g(jbnd,ixkqf(ik,iq0)) / dosef
                      DO iw = 1, nsiw(itemp) ! loop over omega
                         DO iwp = 1, nsiw(itemp) ! loop over omega_prime
                            !
                            ! this step is performed at each iter step only for iw=1 
                            IF ( iw .eq. 1 ) THEN
                               esqrt = 1.d0 / sqrt( wsi(iwp)**2.d0 + ADeltaip(jbnd,ixkqf(ik,iq0),iwp)**2.d0 )
                               wesqrt(jbnd,ixkqf(ik,iq0),iwp) = wsi(iwp) * esqrt 
                               desqrt(jbnd,ixkqf(ik,iq0),iwp) = ADeltaip(jbnd,ixkqf(ik,iq0),iwp) * esqrt 
                            ENDIF
                            IF ( limag_fly ) THEN 
                               CALL lambdar_aniso_ver1( ik, iq, ibnd, jbnd, wsi(iw) - wsi(iwp), lambdam )
                               CALL lambdar_aniso_ver1( ik, iq, ibnd, jbnd, wsi(iw) + wsi(iwp), lambdap )
                            ELSE
                               lambdam = AKeri( ik, iq, ibnd, jbnd, abs(iw-iwp)+1 )
                               lambdap = AKeri( ik, iq, ibnd, jbnd, abs(iw+iwp) )
                            ENDIF
                            kernelm = lambdam - lambdap
                            kernelp = lambdam + lambdap
                            NAZnormi(ibnd,ik,iw) = NAZnormi(ibnd,ik,iw) + weight * kernelm
                            AZnormi(ibnd,ik,iw) = AZnormi(ibnd,ik,iw) + weight * wesqrt(jbnd,ixkqf(ik,iq0),iwp) &
                                                * kernelm
                            ADeltai(ibnd,ik,iw) = ADeltai(ibnd,ik,iw) + weight * desqrt(jbnd,ixkqf(ik,iq0),iwp) &
                                                * ( kernelp - 2.d0 * muc )
                         ENDDO ! iwp
                      ENDDO ! iw
                   ENDIF
                ENDDO ! jbnd
             ENDDO ! iq
          ENDIF
       ENDDO ! ibnd
    ENDDO ! ik
    !
    IF( ALLOCATED(wesqrt) ) DEALLOCATE(wesqrt)
    IF( ALLOCATED(desqrt) ) DEALLOCATE(desqrt)
    !
    ! collect contributions from all pools 
    CALL mp_sum( AZnormi, inter_pool_comm )
    CALL mp_sum( NAZnormi, inter_pool_comm )
    CALL mp_sum( ADeltai, inter_pool_comm )
    CALL mp_barrier(inter_pool_comm)
    !
    IF (mpime .eq. ionode_id) THEN
      IF ( iter .eq. 1 ) THEN
         IF ( .not. ALLOCATED(Deltaold) ) ALLOCATE( Deltaold(nsiw(itemp)) )
         Deltaold(:) = gap0
      ENDIF
      absdelta = 0.d0
      reldelta = 0.d0
      DO iw = 1, nsiw(itemp) ! loop over omega
         DO ik = 1, nkfs
            DO ibnd = 1, nbndfs
               IF ( abs( ekfs(ibnd,ik) - ef0 ) .lt. fsthick ) THEN
                  weight = 0.5d0 * wkfs(ik) * w0g(ibnd,ik) / dosef
                  Znormi(iw) = Znormi(iw) + weight * AZnormi(ibnd,ik,iw)
                  Deltai(iw) = Deltai(iw) + weight * ADeltai(ibnd,ik,iw)
                  NAZnormi(ibnd,ik,iw) = 1.d0 + pi * estemp(itemp) * NAZnormi(ibnd,ik,iw) / wsi(iw)
                  AZnormi(ibnd,ik,iw) = 1.d0 + pi * estemp(itemp) * AZnormi(ibnd,ik,iw) / wsi(iw)
                  ADeltai(ibnd,ik,iw) = pi * estemp(itemp) * ADeltai(ibnd,ik,iw) / AZnormi(ibnd,ik,iw)
               ENDIF
            ENDDO ! ibnd
         ENDDO ! ik
         NZnormi(iw) = 1.d0 + pi * estemp(itemp) * NZnormi(iw) / wsi(iw)
         Znormi(iw) = 1.d0 + pi * estemp(itemp) * Znormi(iw) / wsi(iw)
         Deltai(iw) = pi * estemp(itemp) * Deltai(iw) / Znormi(iw)
         reldelta = reldelta + abs( Deltai(iw) - Deltaold(iw) )
         absdelta = absdelta + abs( Deltai(iw) )
      ENDDO ! iw
      errdelta = reldelta / absdelta
      Deltaold(:) = Deltai(:)
      !
      WRITE(stdout,'(5x,a,i6,a,ES20.10,a,ES20.10,a,ES20.10,a,ES20.10)') 'iter = ', iter, & 
                   '   relerr = ', errdelta, '   abserr = ', reldelta / dble(nsiw(itemp)), &
                   '   Znormi(1) = ', Znormi(1), '   Deltai(1) = ', Deltai(1)
      !
      IF ( errdelta .lt. conv_thr_iaxis) conv = .true.
      IF ( errdelta .lt. conv_thr_iaxis .OR. iter .eq. nsiter ) THEN
         gap(itemp) = Deltai(1)
         gap0 = gap(itemp)
         !
         CALL eliashberg_write_iaxis( itemp )
         !
      ENDIF
      !
      IF ( conv .OR. iter .eq. nsiter ) THEN
         IF( ALLOCATED(Deltaold) ) DEALLOCATE(Deltaold)
         WRITE(stdout,'(5x,a,i6)') 'Convergence was reached in nsiter = ', iter
      ENDIF
      IF ( .not. conv .AND. iter .eq. nsiter ) THEN
         WRITE(stdout,'(a)') ' '
         WRITE(stdout,'(5x,a,i6)') 'Convergence was not reached in nsiter = ', iter
         CALL errore('sum_eliashberg_aniso_iaxis','increase nsiter or reduce conv_thr_iaxis',1)
      ENDIF
    ENDIF
    CALL mp_bcast( Deltai, ionode_id, inter_pool_comm )
    CALL mp_bcast( Znormi, ionode_id, inter_pool_comm )
    CALL mp_bcast( NZnormi, ionode_id, inter_pool_comm )
    CALL mp_bcast( AZnormi, ionode_id, inter_pool_comm )
    CALL mp_bcast( NAZnormi, ionode_id, inter_pool_comm )
    CALL mp_bcast( gap0, ionode_id, inter_pool_comm )
    CALL mp_bcast( gap, ionode_id, inter_pool_comm )
    CALL mp_bcast( Agap, ionode_id, inter_pool_comm )
    CALL mp_bcast( conv, ionode_id, inter_pool_comm )
    CALL mp_barrier(inter_pool_comm)
    !
    IF ( conv .OR. iter .eq. nsiter ) THEN 
       !
       ! remove memory allocated for wesqrt, desqrt, ADeltaip, Deltaold
       imelt = ( 1 + 3 * nbndfs * nkfs ) * nsiw(itemp)
       CALL mem_size_eliashberg( -imelt )
       !
       IF ( .not. limag_fly ) THEN
          !
          IF ( ALLOCATED(AKeri) ) DEALLOCATE(AKeri)
          !
          ! remove memory allocated for AKeri 
          imelt = ( upper_bnd - lower_bnd + 1 ) * maxval(nqfs(:)) * nbndfs**2 * ( 2 * nsiw(itemp) )
          CALL mem_size_eliashberg( -imelt )
          !
       ENDIF
       !
    ENDIF
    !
    RETURN
    !
    END SUBROUTINE sum_eliashberg_aniso_iaxis
    !
    !-----------------------------------------------------------------------
    SUBROUTINE eliashberg_read_aniso_iaxis( itemp )
    !-----------------------------------------------------------------------
    !!  
    !! This routine reads from file the anisotropic Delta and Znorm on the imaginary-axis
    !! 
    !! input
    !!
    !! itemp  - temperature point
    !!
    !---------------------------------------------------------------------- 
    USE kinds,         ONLY : DP
    USE io_epw,        ONLY : iufilgap
    USE io_files,      ONLY : prefix
    USE control_flags, ONLY : iverbosity
    USE epwcom,        ONLY : nstemp, fsthick
    USE eliashbergcom, ONLY : nsiw, estemp, gap0, gap, Agap, wsi, NZnormi, Znormi, Deltai, & 
                              AZnormi, NAZnormi, ADeltai, nkfs, nbndfs, ef0, ekfs, &
                              dosef, wkfs, w0g
    USE constants_epw, ONLY : kelvin2eV, eps6
    USE io_global, ONLY : ionode_id
    USE mp_global, ONLY : inter_pool_comm
    USE mp_world,  ONLY : mpime
    USE mp,        ONLY : mp_bcast, mp_barrier, mp_sum
    ! 
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: itemp
    !
    ! Local variables
    INTEGER :: iw
    !! Counter on frequency
    INTEGER :: ik
    !! Counter on k-poin
    INTEGER :: ibnd
    !! Counter on band
    INTEGER :: imelt
    !! Required allocation of memory
    INTEGER :: ios
    !! Status variables when reading a file
    REAL(DP) :: temp, eband, omega, weight
    CHARACTER (len=256) :: name1, word
    !
    ! get the size of required allocated memory 
    imelt = ( 1 + nbndfs * nkfs ) * nstemp + ( 3 + 3 * nbndfs * nkfs ) * nsiw(itemp)
    CALL mem_size_eliashberg( imelt )
    !
    IF ( .not. ALLOCATED(gap) )      ALLOCATE( gap(nstemp) )
    IF ( .not. ALLOCATED(Agap) )     ALLOCATE( Agap(nbndfs,nkfs,nstemp) )
    IF ( .not. ALLOCATED(Deltai) )   ALLOCATE( Deltai(nsiw(itemp)) )
    IF ( .not. ALLOCATED(Znormi) )   ALLOCATE( Znormi(nsiw(itemp)) )
    IF ( .not. ALLOCATED(NZnormi) )  ALLOCATE( NZnormi(nsiw(itemp)) )
    IF ( .not. ALLOCATED(ADeltai) )  ALLOCATE( ADeltai(nbndfs,nkfs,nsiw(itemp)) )
    IF ( .not. ALLOCATED(AZnormi) )  ALLOCATE( AZnormi(nbndfs,nkfs,nsiw(itemp)) )
    IF ( .not. ALLOCATED(NAZnormi) ) ALLOCATE( NAZnormi(nbndfs,nkfs,nsiw(itemp)) )
    gap(:) = 0.d0
    Agap(:,:,:) = 0.d0
    Deltai(:) = 0.d0
    Znormi(:) = 0.d0
    NZnormi(:) = 0.d0
    ADeltai(:,:,:) = 0.d0
    AZnormi(:,:,:) = 0.d0
    NAZnormi(:,:,:) = 0.d0
    !
    IF (mpime .eq. ionode_id) THEN     
      !   
      temp = estemp(itemp) / kelvin2eV
      ! anisotropic case
      IF ( temp .lt. 10.d0 ) THEN
         WRITE(name1,'(a,a14,f4.2)') TRIM(prefix),'.imag_aniso_00', temp
      ELSEIF ( temp .ge. 10.d0 ) THEN
         WRITE(name1,'(a,a13,f5.2)') TRIM(prefix),'.imag_aniso_0', temp
      ELSEIF ( temp .ge. 100.d0 ) THEN
         WRITE(name1,'(a,a12,f6.2)') TRIM(prefix),'.imag_aniso_', temp
      ENDIF 
      ! 
      OPEN(iufilgap, file=name1, form='formatted', err=100, iostat=ios)
100 CALL errore('eliashberg_read_aniso_iaxis','opening file '//name1,abs(ios))
      READ(iufilgap,'(a)') word
      DO iw = 1, nsiw(itemp) ! loop over omega
         DO ik = 1, nkfs
            DO ibnd = 1, nbndfs
               IF ( abs( ekfs(ibnd,ik) - ef0 ) .lt. fsthick ) THEN
                  READ(iufilgap,'(5ES20.10)') omega, eband, AZnormi(ibnd,ik,iw), ADeltai(ibnd,ik,iw), NAZnormi(ibnd,ik,iw)
                  IF ( iw .eq. 1 ) & 
                     Agap(ibnd,ik,itemp) = ADeltai(ibnd,ik,1)
               ENDIF
            ENDDO ! ibnd
         ENDDO ! ik             
         IF ( abs(wsi(iw)-omega) .gt. eps6 ) &
            CALL errore('eliashberg_read_aniso_iaxis','temperature not the same with the input',1)
      ENDDO ! iw
      CLOSE(iufilgap)
      !
      DO iw = 1, nsiw(itemp) ! loop over omega
        DO ik = 1, nkfs
           DO ibnd = 1, nbndfs
              IF ( abs( ekfs(ibnd,ik) - ef0 ) .lt. fsthick ) THEN
                 weight = 0.5d0 * wkfs(ik) * w0g(ibnd,ik) / dosef
                 Znormi(iw) = Znormi(iw) + weight * AZnormi(ibnd,ik,iw)
                 Deltai(iw) = Deltai(iw) + weight * ADeltai(ibnd,ik,iw)
                 NZnormi(iw) = NZnormi(iw) + weight * NAZnormi(ibnd,ik,iw)
              ENDIF
           ENDDO ! ibnd
        ENDDO ! ik
      ENDDO ! iw
      gap(itemp) = Deltai(1)
      gap0 = gap(itemp)
      !
      CALL gap_FS( itemp )
      !
      IF ( iverbosity .eq. 2 ) &
         CALL free_energy( itemp )
      !
    ENDIF
    CALL mp_bcast( Deltai, ionode_id, inter_pool_comm )
    CALL mp_bcast( Znormi, ionode_id, inter_pool_comm )
    CALL mp_bcast( NZnormi, ionode_id, inter_pool_comm )
    CALL mp_bcast( ADeltai, ionode_id, inter_pool_comm )
    CALL mp_bcast( AZnormi, ionode_id, inter_pool_comm )
    CALL mp_bcast( NAZnormi, ionode_id, inter_pool_comm )
    CALL mp_bcast( gap0, ionode_id, inter_pool_comm )
    CALL mp_bcast( gap, ionode_id, inter_pool_comm )
    CALL mp_bcast( Agap, ionode_id, inter_pool_comm )
    CALL mp_barrier(inter_pool_comm)
    !
    RETURN
    !
    END SUBROUTINE eliashberg_read_aniso_iaxis
    !
    !-----------------------------------------------------------------------
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
    SUBROUTINE eliashberg_write_iaxis( itemp )
    !-----------------------------------------------------------------------
    !!
    !! This routine writes to files results from the solutions of the Eliashberg equations
    !! on the imaginary-axis
    !!
    USE kinds,         ONLY : DP
    USE io_epw,        ONLY : iufilgap
    USE io_files,      ONLY : prefix
    USE control_flags, ONLY : iverbosity
    USE epwcom,        ONLY : fsthick, laniso, liso
    USE eliashbergcom, ONLY : nsiw, estemp, Agap, wsi, & 
                              NAZnormi, AZnormi, ADeltai, NZnormi, Znormi, & 
                              Deltai, nkfs, nbndfs, ef0, ekfs
    USE constants_epw, ONLY : kelvin2eV 
    ! 
    IMPLICIT NONE
    !
    INTEGER, INTENT (in) :: itemp
    !! Counter for temperature
    !
    ! Local variables
    INTEGER  :: iw, ik, ibnd
    REAL(DP) :: temp
    CHARACTER (len=256) :: name1, cname
    !
    temp = estemp(itemp) / kelvin2eV
    !
    cname = 'imag'
    !
    IF ( laniso ) THEN 
       !
       IF ( temp .lt. 10.d0 ) THEN
          WRITE(name1,'(a,a1,a4,a9,f4.2)') TRIM(prefix), '.', cname, '_aniso_00', temp
       ELSEIF ( temp .ge. 10.d0 .AND. temp .lt. 100.d0 ) THEN
          WRITE(name1,'(a,a1,a4,a8,f5.2)') TRIM(prefix), '.', cname, '_aniso_0', temp
       ELSEIF ( temp .ge. 100.d0 ) THEN
          WRITE(name1,'(a,a1,a4,a7,f6.2)') TRIM(prefix), '.', cname, '_aniso_',temp
       ENDIF     
       OPEN(iufilgap, file=name1, form='formatted')
       WRITE(iufilgap,'(5a20)') '#        w [eV]', 'Enk-Ef [eV]', 'Znorm(w)', 'Delta(w) [eV]', 'NZnorm(w)'
       DO iw = 1, nsiw(itemp) ! loop over omega
          DO ik = 1, nkfs
             DO ibnd = 1, nbndfs
                IF ( abs( ekfs(ibnd,ik) - ef0 ) .lt. fsthick ) THEN
                   WRITE(iufilgap,'(5ES20.10)') wsi(iw), ekfs(ibnd,ik)-ef0,&
                         AZnormi(ibnd,ik,iw), ADeltai(ibnd,ik,iw), NAZnormi(ibnd,ik,iw)
                   IF ( iw .eq. 1 ) Agap(ibnd,ik,itemp) = ADeltai(ibnd,ik,iw)
                ENDIF
             ENDDO ! ibnd                   
          ENDDO ! ik
       ENDDO ! iw
       CLOSE(iufilgap)
       !
       CALL gap_distribution_FS ( itemp, cname )
       !
       CALL gap_FS ( itemp )
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
       WRITE(iufilgap,'(4a20)') 'w [eV]', 'Znorm(w)', 'Delta(w) [eV]', 'NZnorm(w)'
       DO iw = 1, nsiw(itemp) ! loop over omega
          WRITE(iufilgap,'(4ES20.10)') wsi(iw), Znormi(iw), Deltai(iw), NZnormi(iw)
       ENDDO
       CLOSE(iufilgap)
    ENDIF 
    !
    RETURN
    !
    END SUBROUTINE eliashberg_write_iaxis
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
    !                                                                            
    !-----------------------------------------------------------------------
    SUBROUTINE kernel_aniso_iaxis( itemp )
    !-----------------------------------------------------------------------
    !!  
    !! Compute kernels K_{+}(ik,iq,ibnd,jbnd;n,n',T) and K_{-}(ik,iq,ibnd,jbnd;n,n',T)
    !! and store them in memory
    !!
    USE kinds,         ONLY : DP
    USE epwcom,        ONLY : fsthick
    USE eliashbergcom, ONLY : nkfs, nbndfs, nsiw, estemp, AKeri, ekfs, ef0, ixkqf, ixqfs, nqfs
    USE constants_epw, ONLY : pi 
    ! 
    IMPLICIT NONE
    !
    INTEGER, INTENT (in) :: itemp
    !! Counter on temperature
    !
    ! Local variables
    INTEGER  :: ik, iq, iq0, ibnd, jbnd, iw, n, lower_bnd, upper_bnd 
    REAL(DP) :: omega, lambda_eph
    !
    CALL fkbounds( nkfs, lower_bnd, upper_bnd )
    IF ( .not. ALLOCATED(AKeri) ) ALLOCATE( AKeri(lower_bnd:upper_bnd,maxval(nqfs(:)),nbndfs,nbndfs,2*nsiw(itemp)) )
    AKeri(:,:,:,:,:) = 0.d0
    !
    ! RM - if lambdar_aniso_ver2 is used then one needs to CALL evaluate_a2fij
    !
    DO ik = lower_bnd, upper_bnd
       DO ibnd = 1, nbndfs
          IF ( abs( ekfs(ibnd,ik) - ef0 ) .lt. fsthick ) THEN
             DO iq = 1, nqfs(ik)
                ! iq0 - index of q-point on the full q-mesh
                iq0 = ixqfs(ik,iq)
                DO jbnd = 1, nbndfs
                   IF ( abs( ekfs(jbnd,ixkqf(ik,iq0)) - ef0 ) .lt. fsthick ) THEN
                      DO iw = 1, 2*nsiw(itemp)
                         n = iw - 1
                         omega = dble(2*n) * pi * estemp(itemp)
                         CALL lambdar_aniso_ver1( ik, iq, ibnd, jbnd, omega, lambda_eph )
                         !CALL lambdar_aniso_ver2( ik, iq, ibnd, jbnd, omega, lambda_eph )
                         AKeri(ik,iq,ibnd,jbnd,iw) = lambda_eph
                      ENDDO ! iw
                   ENDIF
                ENDDO ! jbnd
             ENDDO ! iq
          ENDIF
       ENDDO ! ibnd
    ENDDO ! ik
    !
    RETURN
    !
    END SUBROUTINE kernel_aniso_iaxis
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
    SUBROUTINE lambdar_aniso_ver1( ik, iq, ibnd, jbnd, omega, lambda_eph )
    !-----------------------------------------------------------------------
    !
    ! computes lambda(ik,iq,ibnd,jbnd;n-n')   
    ! reference H. Choi et. al, Physica C 385, 66 (2003)
    !
    ! input
    !
    ! ik - index k-point 
    ! iq - index q-point 
    ! ibnd - index band ibnd at k-point
    ! jbnd - index band jbnd at k+q-point
    ! omega  - frequency 
    !
    ! output
    !
    ! lampda_eph - electron-phonon coupling lambda_ij(k,k+q;n-n')
    !
    !
    USE kinds,         ONLY : DP
    USE phcom,         ONLY : nmodes
    USE elph2,         ONLY : wf
    USE epwcom,        ONLY : eps_acustic
    USE eliashbergcom, ONLY : ixqfs, g2, dosef
    ! 
    IMPLICIT NONE
    !
    INTEGER :: ik, iq, iq0, ibnd, jbnd, imode
    REAL(DP) :: omega, lambda_eph
    !
    ! iq0 - index of q-point on the full q-mesh
    iq0 = ixqfs(ik,iq)
    lambda_eph = 0.d0
    DO imode = 1, nmodes  ! loop over frequency modes
       IF ( wf(imode,iq0) .gt. eps_acustic ) THEN 
          lambda_eph = lambda_eph + g2(ik,iq,ibnd,jbnd,imode) * wf(imode,iq0) & 
                     / ( wf(imode,iq0)**2.d0 + omega**2.d0 )
       ENDIF
    ENDDO 
    lambda_eph = 2.d0 * lambda_eph * dosef
    !
    RETURN
    !
    END SUBROUTINE lambdar_aniso_ver1
    !
    !-----------------------------------------------------------------------
    SUBROUTINE lambdar_aniso_ver2( ik, iq, ibnd, jbnd, omega, lambda_eph )
    !-----------------------------------------------------------------------
    !
    ! computes lambda(ik,iq,ibnd,jbnd;n-n')   
    ! reference H. Choi et. al, Physica C 385, 66 (2003)
    !
    ! input
    !
    ! ik - index k-point 
    ! iq - index q-point 
    ! ibnd - index band ibnd at k-point
    ! jbnd - index band jbnd at k+q-point
    ! omega  - frequency 
    !
    ! output
    !
    ! lampda_eph - electron-phonon coupling lambda_ij(k,k+q;n-n')
    !
    USE kinds,         ONLY : DP
    USE epwcom,        ONLY : nqstep
    USE eliashbergcom, ONLY : a2fij, wsph, dwsph
    !                 
    IMPLICIT NONE        
    !                    
    INTEGER :: ik, iq, ibnd, jbnd, iwph
    REAL(DP) :: omega, lambda_eph 
    !                          
    lambda_eph = 0.d0
    DO iwph = 1, nqstep
       lambda_eph = lambda_eph + wsph(iwph) * a2fij(ik,iq,ibnd,jbnd,iwph) & 
                  / ( wsph(iwph)**2.d0 + omega**2.d0 )
    ENDDO 
    lambda_eph = 2.d0 * lambda_eph * dwsph
    !
    RETURN
    !
    END SUBROUTINE lambdar_aniso_ver2
    !
    !-----------------------------------------------------------------------
    SUBROUTINE kernel_aniso_iaxis_analytic_cont( itemp )
    !-----------------------------------------------------------------------
    !!  
    !! computes kernels K_{+}(w,iw_n,T) and K_{-}(w,iw_n,T)
    !! reference F. Masiglio, M. Schossmann, and J. Carbotte, PRB 37, 4965 (1988)
    !!
    USE kinds,         ONLY : DP
    USE elph2,         ONLY : wqf
    USE epwcom,        ONLY : muc, fsthick
    USE eliashbergcom, ONLY : nsw, nsiw, ws, wsi, ADeltai, nkfs, nbndfs, dosef, ixkqf, ixqfs, nqfs, & 
                              w0g, ekfs, ef0, ADsumi, AZsumi
    USE mp,            ONLY : mp_bcast, mp_barrier, mp_sum
    !      
    IMPLICIT NONE
    !
    INTEGER :: iw, iwp, itemp, ik, iq, iq0, ibnd, jbnd, lower_bnd, upper_bnd, imelt
    REAL(DP) :: esqrt, kernelp, kernelm, weight
    REAL(DP), ALLOCATABLE :: wesqrt(:,:,:), desqrt(:,:,:)
    COMPLEX(DP) :: lambda_eph
    !
    CALL fkbounds( nkfs, lower_bnd, upper_bnd )
    !
    ! get memory size required for wesqrt, desqrt, ADsumi, AZsumi
    imelt = 2 * nbndfs * nkfs * nsiw(itemp) + 2 * ( upper_bnd - lower_bnd + 1 ) * nbndfs * nsw 
    CALL mem_size_eliashberg( imelt )
    !
    IF ( .not. ALLOCATED(wesqrt) ) ALLOCATE( wesqrt(nbndfs,nkfs,nsiw(itemp)) )
    IF ( .not. ALLOCATED(desqrt) ) ALLOCATE( desqrt(nbndfs,nkfs,nsiw(itemp)) )
    !
    DO ik = lower_bnd, upper_bnd
       IF ( .not. ALLOCATED(ADsumi) ) ALLOCATE( ADsumi(nbndfs,lower_bnd:upper_bnd,nsw) )
       IF ( .not. ALLOCATED(AZsumi) ) ALLOCATE( AZsumi(nbndfs,lower_bnd:upper_bnd,nsw) )
    ENDDO
    ADsumi(:,:,:) = 0.d0
    AZsumi(:,:,:) = 0.d0
    !
    ! RM - if lambdai_aniso_ver2 is used then one needs to CALL evaluate_a2fij
    !
    DO ik = lower_bnd, upper_bnd
       DO ibnd = 1, nbndfs
          IF ( abs( ekfs(ibnd,ik) - ef0 ) .lt. fsthick ) THEN
             DO iq = 1, nqfs(ik)
                ! iq0 - index of q-point on the full q-mesh
                iq0 = ixqfs(ik,iq)
                DO jbnd = 1, nbndfs
                   IF ( abs( ekfs(jbnd,ixkqf(ik,iq0)) - ef0 ) .lt. fsthick ) THEN
                      weight = wqf(iq) * w0g(jbnd,ixkqf(ik,iq0)) / dosef
                      DO iw = 1, nsw ! loop over omega 
                         DO iwp = 1, nsiw(itemp) ! loop over iw_n
                            CALL lambdai_aniso_ver1( ik, iq, ibnd, jbnd, ws(iw), wsi(iwp), lambda_eph )
                            !CALL lambdai_aniso_ver2( ik, iq, ibnd, jbnd, ws(iw), wsi(iwp), lambda_eph )
                            kernelp = 2.d0 * real(lambda_eph)
                            kernelm = 2.d0 * aimag(lambda_eph)
                            IF ( iw .eq. 1 ) THEN
                               esqrt = 1.d0 / sqrt( wsi(iwp)**2.d0 + ADeltai(jbnd,ixkqf(ik,iq0),iwp)**2.d0 )
                               wesqrt(jbnd,ixkqf(ik,iq0),iwp) =  wsi(iwp) * esqrt
                               desqrt(jbnd,ixkqf(ik,iq0),iwp) =  ADeltai(jbnd,ixkqf(ik,iq0),iwp) * esqrt
                            ENDIF
                            AZsumi(ibnd,ik,iw) = AZsumi(ibnd,ik,iw) & 
                                               + weight * wesqrt(jbnd,ixkqf(ik,iq0),iwp) * kernelm 
                            ADsumi(ibnd,ik,iw) = ADsumi(ibnd,ik,iw) & 
                                               + weight * desqrt(jbnd,ixkqf(ik,iq0),iwp) * ( kernelp - 2.d0 * muc ) 
                         ENDDO ! iwp
                      ENDDO ! iw
                   ENDIF
                ENDDO ! jbnd
             ENDDO ! iq
          ENDIF
       ENDDO ! ibnd
    ENDDO ! ik
    !
    IF( ALLOCATED(wesqrt) ) DEALLOCATE(wesqrt)
    IF( ALLOCATED(desqrt) ) DEALLOCATE(desqrt)
    !  
    ! remove memory allocated for wesqrt, desqrt
    imelt = 2 * nbndfs * nkfs * nsiw(itemp) 
    CALL mem_size_eliashberg ( -imelt )
    !
    RETURN
    !
    END SUBROUTINE kernel_aniso_iaxis_analytic_cont
    !                                                
    !-----------------------------------------------------------------------
    SUBROUTINE lambdai_aniso_ver1( ik, iq, ibnd, jbnd, omega, omegap, lambda_eph )
    !-----------------------------------------------------------------------
    !!
    !! computes lambda(w-iw_n)   
    !! reference F. Masiglio, M. Schossmann, and J. Carbotte, PRB 37, 4965 (1988)
    !!
    !! input
    !!  
    !! ik - index k-point
    !! iq - index q-point 
    !! ibnd - index band ibnd at k-point
    !! jbnd - index band jbnd at k+q-point
    !! iw     - index frequency w on the real-axis
    !! iwp    - index frequency iw_n on the imaginary-axis
    !! omega  - frequency w at point iw
    !! omegap - frequency w_n at point iwp
    !!     
    !! output 
    !!        
    !! lampda_eph - electron-phonon coupling lambda_ij(k,k+q;w-iw_n)
    !!        
    USE kinds,         ONLY : DP
    USE phcom,         ONLY : nmodes
    USE elph2,         ONLY : wf
    USE epwcom,        ONLY : eps_acustic
    USE eliashbergcom, ONLY : ixqfs, g2, dosef
    USE constants_epw, ONLY : ci
    !     
    IMPLICIT NONE
    !  
    INTEGER :: ik, iq, iq0, ibnd, jbnd, imode
    REAL(DP) :: omega, omegap
    COMPLEX(DP) :: lambda_eph
    !
    ! iq0 - index of q-point on the full q-mesh
    iq0 = ixqfs(ik,iq)
    lambda_eph = (0.d0,0.d0)
    DO imode = 1, nmodes  ! loop over frequency modes
       IF ( wf(imode,iq0) .gt. eps_acustic ) THEN 
          lambda_eph = lambda_eph +  g2(ik,iq,ibnd,jbnd,imode) * wf(imode,iq0) & 
                     / ( wf(imode,iq0)**2.d0 - (omega - ci*omegap)**2.d0 )
       ENDIF
    ENDDO ! iwph
    lambda_eph = 2.d0 * lambda_eph * dosef
    !
    RETURN
    !
    END SUBROUTINE lambdai_aniso_ver1
    !
    !-----------------------------------------------------------------------               
    SUBROUTINE lambdai_aniso_ver2( ik, iq, ibnd, jbnd, omega, omegap, lambda_eph )
    !-----------------------------------------------------------------------
    !!
    !! computes lambda(w-iw_n)   
    !! reference F. Masiglio, M. Schossmann, and J. Carbotte, PRB 37, 4965 (1988)
    !!
    !! input
    !!  
    !! ik - index k-point
    !! iq - index q-point 
    !! ibnd - index band ibnd at k-point
    !! jbnd - index band jbnd at k+q-point
    !! iw     - index frequency w on the real-axis
    !! iwp    - index frequency iw_n on the imaginary-axis
    !! omega  - frequency w at point iw
    !! omegap - frequency w_n at point iwp
    !!     
    !! output 
    !!        
    !! lampda_eph - electron-phonon coupling lambda_ij(k,k+q;w-iw_n)
    !!        
    USE kinds,         ONLY : DP
    USE epwcom,        ONLY : nqstep
    USE eliashbergcom, ONLY : a2fij, dwsph, wsph
    USE constants_epw, ONLY : ci
    !     
    IMPLICIT NONE
    !  
    INTEGER :: ik, iq, ibnd, jbnd, iwph
    REAL(DP) :: omega, omegap
    COMPLEX(DP) :: lambda_eph
    !
    lambda_eph = (0.d0,0.d0)
    DO iwph = 1, nqstep
       lambda_eph = lambda_eph + wsph(iwph) * a2fij(ik,iq,ibnd,jbnd,iwph) &
                  / ( wsph(iwph)**2.d0 - (omega - ci*omegap)**2.d0 )
    ENDDO ! iwph
    lambda_eph = 2.d0 * lambda_eph * dwsph
    !
    RETURN
    !
    END SUBROUTINE lambdai_aniso_ver2
    !
    !-----------------------------------------------------------------------
    SUBROUTINE evaluate_a2fij
    !-----------------------------------------------------------------------
    !!
    !! computes the anisotropic spectral function a2F(k,k',w) 
    !!
    USE kinds,         ONLY : DP
    USE phcom,         ONLY : nmodes
    USE elph2,         ONLY : wf
    USE epwcom,        ONLY : fsthick, eps_acustic, nqstep, degaussq
    USE eliashbergcom, ONLY : nkfs, nbndfs, g2, a2fij, ixkqf, ixqfs, nqfs, ekfs, ef0, & 
                              dosef, wsph
    ! 
    IMPLICIT NONE
    !
    INTEGER :: ik, iq, iq0, iwph, ibnd, jbnd, imode, lower_bnd, upper_bnd
    REAL(DP) :: weight
    REAL(DP), EXTERNAL :: w0gauss
    !
    CALL fkbounds( nkfs, lower_bnd, upper_bnd )
    IF ( .not. ALLOCATED(a2fij) ) ALLOCATE(a2fij(lower_bnd:upper_bnd,maxval(nqfs(:)),nbndfs,nbndfs,nqstep))
    a2fij(:,:,:,:,:) = 0.d0
    !
    DO ik = lower_bnd, upper_bnd 
       DO ibnd = 1, nbndfs
          IF ( abs( ekfs(ibnd,ik) - ef0 ) .lt. fsthick ) THEN
             DO iq = 1, nqfs(ik)
                ! iq0 - index of q-point on the full q-mesh
                iq0 = ixqfs(ik,iq)
                DO jbnd = 1, nbndfs
                   IF ( abs( ekfs(jbnd,ixkqf(ik,iq0)) - ef0 ) .lt. fsthick ) THEN
                      DO imode = 1, nmodes
                         IF ( wf(imode,iq0) .gt. eps_acustic ) THEN 
                            DO iwph = 1, nqstep
                               weight  = w0gauss( ( wsph(iwph) - wf(imode,iq0) ) / degaussq, 0 ) / degaussq
                               a2fij(ik,iq,ibnd,jbnd,iwph) = a2fij(ik,iq,ibnd,jbnd,iwph) &
                                                           + weight * dosef * g2(ik,iq,ibnd,jbnd,imode)
                            ENDDO ! iwph
                         ENDIF ! wf
                      ENDDO ! imode
                   ENDIF
                ENDDO ! jbnd
             ENDDO ! iq
          ENDIF
       ENDDO ! ibnd
    ENDDO ! ik
    !
    RETURN
    !
    END SUBROUTINE evaluate_a2fij
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
    !   
    !-----------------------------------------------------------------------
    !                                                                            
    !-----------------------------------------------------------------------
    SUBROUTINE kernel_iso_iaxis( itemp )
    !-----------------------------------------------------------------------
    !  
    ! computes kernels K_{+}(n,n',T) and K_{-}(n,n'T)
    ! reference W. E. Pickett, PRB 26, 1186 (1982)
    !
    !
    USE kinds, ONLY : DP
    USE constants_epw, ONLY : pi
    USE eliashbergcom, ONLY : nsiw, estemp, Keri
    ! 
    IMPLICIT NONE
    !
    INTEGER  :: iw, itemp, n
    REAL(DP) :: omega, lambda_eph
    !
    IF ( .not. ALLOCATED(Keri) ) ALLOCATE( Keri(2*nsiw(itemp)) )
    Keri(:) = 0.d0
    !
    DO iw = 1, 2*nsiw(itemp)
       n = iw - 1
       omega = dble(2*n) * pi * estemp(itemp)
       CALL lambdar_iso( omega, lambda_eph )
       Keri(iw) = lambda_eph
    ENDDO 
    !
    RETURN
    !
    END SUBROUTINE kernel_iso_iaxis                                                       
    !
    !-----------------------------------------------------------------------
    SUBROUTINE lambdar_iso( omega, lambda_eph )
    !-----------------------------------------------------------------------
    !
    ! computes lambda(n-n')   
    ! reference W. E. Pickett, PRB 26, 1186 (1982)
    !
    ! input
    !
    ! omega  - frequency 
    !
    ! output
    !
    ! lampda_eph - electron-phonon coupling lambda(n-n')
    !
    !
    USE kinds, ONLY : DP
    USE epwcom, ONLY : nqstep
    USE eliashbergcom, ONLY : a2f_iso, wsph, dwsph
    ! 
    IMPLICIT NONE
    !
    INTEGER :: iwph
    REAL(DP) :: omega, lambda_eph
    !
    lambda_eph = 0.d0
    DO iwph = 1, nqstep  ! loop over Omega (integration variable)
       lambda_eph = lambda_eph + wsph(iwph) * a2f_iso(iwph) & 
                  / ( wsph(iwph)**2.d0 + omega**2.d0 )
    ENDDO ! iwph
    lambda_eph = 2.d0 * lambda_eph * dwsph 
    !
    RETURN
    !
    END SUBROUTINE lambdar_iso
  
    !-----------------------------------------------------------------------
    SUBROUTINE kernel_iso_iaxis_analytic_cont( itemp )
    !-----------------------------------------------------------------------
    !  
    ! computes kernels K_{+}(w,iw_n,T) and K_{-}(w,iw_n,T)
    ! reference F. Masiglio, M. Schossmann, and J. Carbotte, PRB 37, 4965 (1988)
    !
    !
    USE kinds,         ONLY : DP
    USE epwcom,        ONLY : muc
    USE eliashbergcom, ONLY : nsw, nsiw, ws, wsi, Deltai, Dsumi, Zsumi
    ! 
    IMPLICIT NONE
    !
    INTEGER :: iw, iwp, itemp
    REAL(DP) :: esqrt, kernelp, kernelm
    REAL(DP), ALLOCATABLE :: wesqrt(:), desqrt(:)
    COMPLEX(DP) :: lambda_eph
    !
    IF ( .not. ALLOCATED(wesqrt) ) ALLOCATE( wesqrt(nsiw(itemp)) )
    IF ( .not. ALLOCATED(desqrt) ) ALLOCATE( desqrt(nsiw(itemp)) )
    IF ( .not. ALLOCATED(Dsumi) )  ALLOCATE( Dsumi(nsw) )
    IF ( .not. ALLOCATED(Zsumi) )  ALLOCATE( Zsumi(nsw) )
    Dsumi(:) = 0.d0
    Zsumi(:) = 0.d0
    !
    DO iw = 1, nsw ! loop over omega
       DO iwp = 1, nsiw(itemp) ! loop over iw_n
          CALL lambdai_iso( ws(iw), wsi(iwp), lambda_eph )
          kernelp = 2.d0 * real(lambda_eph)
          kernelm = 2.d0 * aimag(lambda_eph) 
          IF ( iw .eq. 1 ) THEN
             esqrt = 1.d0 / sqrt( wsi(iwp)**2.d0 + Deltai(iwp)**2.d0 )
             wesqrt(iwp) =  wsi(iwp) * esqrt
             desqrt(iwp) =  Deltai(iwp) * esqrt
          ENDIF
          Zsumi(iw) = Zsumi(iw) + kernelm * wesqrt(iwp)
          Dsumi(iw) = Dsumi(iw) + ( kernelp - 2.d0 * muc ) * desqrt(iwp)
       ENDDO
    ENDDO
    !
    IF( ALLOCATED(wesqrt) ) DEALLOCATE (wesqrt)
    IF( ALLOCATED(desqrt) ) DEALLOCATE (desqrt)
    !   
    RETURN
    !
    END SUBROUTINE kernel_iso_iaxis_analytic_cont      
    !                                                
    !-----------------------------------------------------------------------
    SUBROUTINE lambdai_iso( omega, omegap, lambda_eph )
    !-----------------------------------------------------------------------
    !
    ! computes lambda(w-iw_n)   
    ! reference F. Masiglio, M. Schossmann, and J. Carbotte, PRB 37, 4965 (1988)
    !
    ! input
    !
    ! iw     - index frequency w on the real-axis
    ! iwp    - index frequency iw_n on the imaginary-axis
    ! omega  - frequency w at point iw
    ! omegap - frequency w_n at point iwp
    !
    ! output
    !
    ! lampda_eph - electron-phonon coupling lambda(w-iw_n)
    !
    !
    USE kinds, ONLY : DP
    USE epwcom,        ONLY : nqstep
    USE eliashbergcom, ONLY : a2f_iso, wsph, dwsph
    USE constants_epw, ONLY : ci
    ! 
    IMPLICIT NONE
    !
    INTEGER :: iwph
    REAL(DP) :: omega, omegap
    COMPLEX(DP) :: lambda_eph
    !
    lambda_eph = (0.d0,0.d0)
    DO iwph = 1, nqstep  ! loop over Omega (integration variable)
       lambda_eph = lambda_eph & 
                  + wsph(iwph) * a2f_iso(iwph) / ( wsph(iwph)**2.d0 - (omega - ci*omegap)**2.d0 )
    ENDDO ! iwph
    lambda_eph = lambda_eph * 2.d0 * dwsph 
    !
    RETURN
    !
    END SUBROUTINE lambdai_iso
    !                                        
    !-----------------------------------------------------------------------  
    !                                                                            
    !-----------------------------------------------------------------------
    SUBROUTINE kernel_raxis( iw, iwp, itemp, kernelp, kernelm )
    !-----------------------------------------------------------------------
    !
    ! computes kernels K_{+}(w,w',T) and K_{-}(w,w'T)  
    ! reference M. J. Holcomb, PRB 54, 6648 (1996)   
    !
    ! input
    !
    ! iw     - index frequency w : ws(iw)
    ! iwp    - index frequency w' : ws(iwp)
    ! itemp  - index temperature
    !
    ! output
    !
    ! kernelp - phonon kernel K_{+}(w,w',T)
    ! kernelm - phonon kernel K_{-}(w,w',T)
    !
    USE kinds,         ONLY : DP
    USE constants_epw, ONLY : pi, ci, eps6
    USE epwcom,        ONLY : nqstep
    USE eliashbergcom, ONLY : a2f_iso, bewph, wsph, dwsph, ws, fdwp, estemp
    ! 
    IMPLICIT NONE
    !
    INTEGER :: iw, iwp, iwph, itemp, ngaussw0
    REAL(DP) :: degaussw0, f1, f2, f3, f4
    COMPLEX(DP) :: e1, e2, e3, e4, kernelp, kernelm
    REAL(DP), EXTERNAL :: wgauss, w0gauss
    !
    degaussw0 = 1.d0 * dwsph
    ngaussw0 = 0
    !
    f1 = 0.d0
    f2 = 0.d0
    f3 = 0.d0
    f4 = 0.d0
    kernelp = (0.d0, 0.d0)
    kernelm = (0.d0, 0.d0)
    e1 = (0.d0, 0.d0)
    e2 = (0.d0, 0.d0)
    e3 = (0.d0, 0.d0)
    e4 = (0.d0, 0.d0)
    !
    IF ( .not. ALLOCATED(bewph) ) ALLOCATE( bewph(nqstep) )
    ! Bose-Einstein distribution
    DO iwph = 1, nqstep  ! loop over Omega (integration variable)
       IF ( iw .eq. 1 .AND. iwp .eq. 1 ) THEN
          IF ( ABS(estemp(itemp)) <  eps6 ) THEN
             bewph(iwph)  = 0.d0
          ELSE
             bewph(iwph) = wgauss( -wsph(iwph) / estemp(itemp), -99 )
             bewph(iwph) = bewph(iwph) / (1.d0 - 2.d0 * bewph(iwph))
          ENDIF
       ENDIF
       !
       ! a small complex number is added to denominator to move the pole away from the real-axis
       !
       ! in order to reduce the numerical noise at very small frequencies coming from the complex number 
       ! added in the denominator, the contribution of the imaginary part is reestimated using 
       ! delta function (RM notes) 
       !
       ! subtract the imaginary part coming from e1 to e4 and add instead the imaginary part 
       ! coming from f1 to f4
       !
       e1 = 1.d0 / ( wsph(iwph) + ws(iwp) + ws(iw) + ci*degaussw0 ) 
       e2 = 1.d0 / ( wsph(iwph) + ws(iwp) - ws(iw) - ci*degaussw0 ) 
       e3 = 1.d0 / ( wsph(iwph) - ws(iwp) + ws(iw) + ci*degaussw0 ) 
       e4 = 1.d0 / ( wsph(iwph) - ws(iwp) - ws(iw) - ci*degaussw0 ) 
       !
       ! estimate of the imaginary part using delta function
       f1 = w0gauss( ( wsph(iwph) + ws(iwp) + ws(iw) )/degaussw0, ngaussw0 ) / degaussw0
       f2 = w0gauss( ( wsph(iwph) + ws(iwp) - ws(iw) )/degaussw0, ngaussw0 ) / degaussw0
       f3 = w0gauss( ( wsph(iwph) - ws(iwp) + ws(iw) )/degaussw0, ngaussw0 ) / degaussw0
       f4 = w0gauss( ( wsph(iwph) - ws(iwp) - ws(iw) )/degaussw0, ngaussw0 ) / degaussw0
       !
       kernelp = kernelp + a2f_iso(iwph) &
               * (  ( 1.d0 - fdwp(iwp) + bewph(iwph) ) * ( e1 - ci*aimag(e1) - ci*pi*f1 + e2 - ci*aimag(e2) + ci*pi*f2 ) & 
                  - (        fdwp(iwp) + bewph(iwph) ) * ( e3 - ci*aimag(e3) - ci*pi*f3 + e4 - ci*aimag(e4) + ci*pi*f4 ) )
       kernelm = kernelm + a2f_iso(iwph) &
               * (  ( 1.d0 - fdwp(iwp) + bewph(iwph) ) * ( e1 - ci*aimag(e1) - ci*pi*f1 - e2 + ci*aimag(e2) - ci*pi*f2 ) &
                  + (        fdwp(iwp) + bewph(iwph) ) * ( e3 - ci*aimag(e3) - ci*pi*f3 - e4 + ci*aimag(e4) - ci*pi*f4 ) )
    ENDDO ! iwph
    kernelp = kernelp * dwsph 
    kernelm = kernelm * dwsph 
    !
    RETURN
    !
    END SUBROUTINE kernel_raxis
    !
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
    USE, INTRINSIC :: IEEE_ARITHMETIC
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
      IF ( IEEE_IS_NAN(ar) .or. IEEE_IS_NAN(ai) ) THEN
  !     write(6,*) (z(i),i=1,N)
  !     write(6,*) (u(i),i=1,N)
  !     write(6,*) (a(i),i=1,N)
        write(6,*) 'one or more coefficients are NaN'
  !     call errore('pade_coeff','one or more coefficients are NaN',1)
      ENDIF
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
    !-----------------------------------------------------------
  
    !-----------------------------------------------------------------------
    SUBROUTINE gap_distribution_FS ( itemp, cname )
    !-----------------------------------------------------------------------
    !
    ! This routine writes to files the distribution of the superconducting 
    ! gap on the Fermi surface
    !
    USE kinds,         ONLY : DP
    USE io_epw,        ONLY : iufilgap
    USE io_files,      ONLY : prefix
    USE epwcom,        ONLY : fsthick
    USE eliashbergcom, ONLY : estemp, Agap, nkfs, nbndfs, ef0, ekfs
    USE constants_epw, ONLY : kelvin2eV
    !
    IMPLICIT NONE
    !
    INTEGER  :: ik, ibnd, ibin, nbin, itemp
    REAL(DP) :: weight, temp, delta_max, dbin, sigma
    REAL(DP), ALLOCATABLE :: delta_k_bin(:)
    REAL(DP), EXTERNAL :: w0gauss
    CHARACTER (len=256) :: name1, cname
    !
    temp = estemp(itemp) / kelvin2eV
    !
    delta_max = 1.25d0 * maxval(Agap(:,:,itemp))
    nbin = int(delta_max/(0.005d0/1000.d0))
    dbin = delta_max / dble(nbin)
    IF ( .not. ALLOCATED(delta_k_bin) ) ALLOCATE( delta_k_bin(nbin) )
    delta_k_bin(:) = 0.d0
    !
    DO ik = 1, nkfs
       DO ibnd = 1, nbndfs
          IF ( abs( ekfs(ibnd,ik) - ef0 ) .lt. fsthick ) THEN
             DO ibin = 1, nbin
                sigma = 1.d0 * dbin
                weight = w0gauss( ( Agap(ibnd,ik,itemp) - dble(ibin) * dbin) / sigma, 0 ) / sigma
                delta_k_bin(ibin) = delta_k_bin(ibin) + weight
             ENDDO
          ENDIF
       ENDDO
    ENDDO
    !
    IF ( temp .lt. 10.d0 ) THEN
       WRITE(name1,'(a,a1,a4,a14,f4.2)') TRIM(prefix), '.', cname, '_aniso_gap0_00', temp
    ELSEIF ( temp .ge. 10.d0 .AND. temp .lt. 100.d0  ) THEN
       WRITE(name1,'(a,a1,a4,a13,f5.2)') TRIM(prefix), '.', cname, '_aniso_gap0_0', temp
    ELSEIF ( temp .ge. 100.d0 ) THEN
       WRITE(name1,'(a,a1,a4,a12,f6.2)') TRIM(prefix), '.', cname, '_aniso_gap0_', temp
    ENDIF
    !
    OPEN(iufilgap, file=name1, form='formatted')
    DO ibin = 1, nbin
       WRITE(iufilgap,'(2ES20.10)') temp + delta_k_bin(ibin)/maxval(delta_k_bin(:)), dbin*dble(ibin)
    ENDDO
    CLOSE(iufilgap)
    !
    IF ( ALLOCATED(delta_k_bin) ) DEALLOCATE(delta_k_bin)
    !
    RETURN
    !
    END SUBROUTINE gap_distribution_FS
    !
    !-----------------------------------------------------------------------
    SUBROUTINE gap_FS( itemp )
    !-----------------------------------------------------------------------
    !
    ! This routine writes to files the superconducting gap on the Fermi surface
    !
    USE kinds,         ONLY : DP
    USE io_epw,        ONLY : iufilgapFS
    USE io_files,      ONLY : prefix
    USE cell_base,     ONLY : bg
    USE control_flags, ONLY : iverbosity
    USE epwcom,        ONLY : fsthick, nkf1, nkf2, nkf3
    USE eliashbergcom, ONLY : estemp, Agap, nkfs, nbndfs, ef0, ekfs, ixkff
    USE constants_epw, ONLY : kelvin2eV
    !
    IMPLICIT NONE
    !
    INTEGER  :: i, j, k, itemp, ik, ibnd
    REAL(DP) :: temp, x1, x2, x3
    REAL(DP), ALLOCATABLE :: Agap_tmp(:,:)
    CHARACTER (len=256) :: name1, cname
    !
    temp = estemp(itemp) / kelvin2eV
    !
    cname = 'imag'
    !
    ! RM - If the k-point is outside the Fermi shell,
    ! ixkff(ik)=0 and Agap_tmp(:,0) = 0.0
    !
    IF ( .not. ALLOCATED(Agap_tmp) ) ALLOCATE(Agap_tmp(nbndfs,0:nkfs))
    Agap_tmp(:,1:nkfs) = Agap(:,1:nkfs,itemp)
    Agap_tmp(:,0) = 0.0d0
    !
    ! SP & RM: .cube file for VESTA plotting (only if iverbosity = 2)
    !
    IF ( iverbosity .eq. 2 ) THEN
      !
      DO ibnd = 1, nbndfs
        !
        IF ( ibnd < 10 ) THEN
          ! We make the assumption that there are no superconductor with Tc higher than 999 K.  
          IF ( temp < 10.d0 ) THEN
             WRITE(name1,'(a,a1,a4,a14,f4.2,a1,i1,a5)')TRIM(prefix), '.', cname, '_aniso_gap0_00', temp, '_', ibnd, '.cube'
          ELSEIF ( temp < 100.d0 ) THEN
             WRITE(name1,'(a,a1,a4,a13,f5.2,a1,i1,a5)')TRIM(prefix), '.', cname, '_aniso_gap0_0', temp, '_', ibnd, '.cube'
          ELSEIF ( temp < 1000.d0 ) THEN
             WRITE(name1,'(a,a1,a4,a12,f6.2,a1,i1,a5)')TRIM(prefix), '.', cname, '_aniso_gap0_', temp, '_', ibnd, '.cube'
          ENDIF
        ELSEIF ( ibnd < 100 ) THEN
          IF ( temp < 10.d0 ) THEN
             WRITE(name1,'(a,a1,a4,a14,f4.2,a1,i2,a5)')TRIM(prefix), '.', cname, '_aniso_gap0_00', temp, '_', ibnd, '.cube'
          ELSEIF ( temp < 100.d0 .and. temp > 9.9999d0 ) THEN
             WRITE(name1,'(a,a1,a4,a13,f5.2,a1,i2,a5)')TRIM(prefix), '.', cname, '_aniso_gap0_0', temp, '_', ibnd, '.cube'
          ELSEIF ( temp < 1000.d0 .and. temp > 99.9999d0 ) THEN
             WRITE(name1,'(a,a1,a4,a12,f6.2,a1,i2,a5)')TRIM(prefix), '.', cname, '_aniso_gap0_', temp, '_', ibnd, '.cube'
          ENDIF
        ELSEIF ( ibnd < 1000 ) THEN
          IF ( temp < 10.d0 ) THEN
             WRITE(name1,'(a,a1,a4,a14,f4.2,a1,i3,a5)')TRIM(prefix), '.', cname, '_aniso_gap0_00', temp, '_', ibnd, '.cube'
          ELSEIF ( temp < 100.d0 .and. temp > 9.9999d0  ) THEN
             WRITE(name1,'(a,a1,a4,a13,f5.2,a1,i3,a5)')TRIM(prefix), '.', cname, '_aniso_gap0_0', temp, '_', ibnd, '.cube'
          ELSEIF ( temp < 1000.d0 .and. temp > 99.9999d0 ) THEN
             WRITE(name1,'(a,a1,a4,a12,f6.2,a1,i3,a5)')TRIM(prefix), '.', cname, '_aniso_gap0_', temp, '_', ibnd, '.cube'
          ENDIF
        ELSE
          CALL errore( 'eliashberg_write', ' Too many bands ',1)
        ENDIF
        !
        OPEN(iufilgapFS, file=name1, form='formatted')
        WRITE(iufilgapFS,*) 'Cubfile created from EPW calculation'
        WRITE(iufilgapFS,*) 'gap'
        WRITE(iufilgapFS,'(i5,3f12.6)') 1, 0.0d0, 0.0d0, 0.0d0
        WRITE(iufilgapFS,'(i5,3f12.6)') nkf1, (bg(i,1)/dble(nkf1),i=1,3)
        WRITE(iufilgapFS,'(i5,3f12.6)') nkf2, (bg(i,2)/dble(nkf2),i=1,3)
        WRITE(iufilgapFS,'(i5,3f12.6)') nkf3, (bg(i,3)/dble(nkf3),i=1,3)
        WRITE(iufilgapFS,'(i5,4f12.6)') 1, 1.0d0, 0.0d0, 0.0d0, 0.0d0
        WRITE(iufilgapFS,'(6f12.6)') ( Agap_tmp(ibnd,ixkff(ik)),ik=1,nkf1*nkf2*nkf3 )
        CLOSE(iufilgapFS)
      ENDDO
      !
    ENDIF
    !
    ! SP & RM : Write on file the superconducting gap close to the Fermi surface along with
    !     Cartesian coordinate, band index, energy distance from Fermi level and gap value.
    !
    IF ( temp .lt. 10.d0 ) THEN
       WRITE(name1,'(a,a1,a4,a16,f4.2)') TRIM(prefix), '.', cname, '_aniso_gap_FS_00', temp
    ELSEIF ( temp .ge. 10.d0 .AND. temp .lt. 100.d0) THEN
       WRITE(name1,'(a,a1,a4,a15,f5.2)') TRIM(prefix), '.', cname, '_aniso_gap_FS_0', temp
    ELSEIF ( temp .ge. 100.d0 ) THEN
       WRITE(name1,'(a,a1,a4,a14,f6.2)') TRIM(prefix), '.', cname, '_aniso_gap_FS_', temp
    ENDIF
    OPEN(iufilgapFS, file=name1, form='formatted')
    WRITE(iufilgapFS,'(a78)') '#               k-point                  Band Enk-Ef [eV]        Delta(0) [eV]'
    DO i = 1, nkf1
      DO j = 1, nkf2
        DO k = 1, nkf3
          ik = k + (j-1)*nkf3 + (i-1)*nkf2*nkf3
          !IF ( ixkff(ik) .gt. 0 ) THEN
            DO ibnd = 1, nbndfs
              ! RM: Everything is in eV here.
              ! SP: Here take a 0.2 eV interval around the FS.
              IF ( abs( ekfs(ibnd,ixkff(ik)) - ef0 ) .lt. fsthick ) THEN
              !IF ( abs( ekfs(ibnd,ixkff(ik)) - ef0 ) .lt. 0.2 ) THEN
                 x1 = bg(1,1)*(i-1)/nkf1+bg(1,2)*(j-1)/nkf2+bg(1,3)*(k-1)/nkf3
                 x2 = bg(2,1)*(i-1)/nkf1+bg(2,2)*(j-1)/nkf2+bg(2,3)*(k-1)/nkf3
                 x3 = bg(3,1)*(i-1)/nkf1+bg(3,2)*(j-1)/nkf2+bg(3,3)*(k-1)/nkf3
                 WRITE(iufilgapFS,'(3f12.6,i8,f12.6,f24.15)') x1, x2, x3, ibnd, &
                       ekfs(ibnd,ixkff(ik))-ef0, Agap_tmp(ibnd,ixkff(ik))
              ENDIF
            ENDDO ! ibnd
          !ENDIF
        ENDDO  ! k
      ENDDO ! j
    ENDDO ! i
    CLOSE(iufilgapFS)
    !
    IF ( ALLOCATED(Agap_tmp) ) DEALLOCATE(Agap_tmp)
    !
    RETURN
    !
    END SUBROUTINE gap_FS
    !
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE read_a2f
    !-----------------------------------------------------------------------
    !!
    !! Read the eliashberg spectral function from fila2f
    !!
    USE io_global,     ONLY : stdout
    USE epwcom,        ONLY : nqstep, fila2f
    USE eliashbergcom, ONLY : wsphmax, wsph, a2f_iso, memlt_pool
    USE mp_global,     ONLY : npool
    USE io_epw,        ONLY : iua2ffil 
    USE io_global, ONLY : ionode_id
    USE mp_global, ONLY : inter_pool_comm
    USE mp,        ONLY : mp_bcast, mp_barrier, mp_sum
    USE mp_world,  ONLY : mpime
    ! 
    IMPLICIT NONE
    !
    INTEGER :: iwph
    !! Counter for the number of freq
    INTEGER :: ios
    !! Status when opening a2F file
    !
    IF ( .not. ALLOCATED(a2f_iso) ) ALLOCATE(a2f_iso(nqstep))
    IF ( .not. ALLOCATED(wsph) ) ALLOCATE(wsph(nqstep)) 
    a2f_iso(:) = 0.d0
    wsph(:) = 0.d0
    !
    IF ( mpime .eq. ionode_id ) THEN
      OPEN(iua2ffil, file=fila2f, status='unknown', err=100, iostat=ios)
100   CALL errore('read_a2f','opening file'//fila2f,abs(ios))
    !
      DO iwph = 1, nqstep
         READ(iua2ffil,*) wsph(iwph), a2f_iso(iwph) ! freq from meV to eV
         wsph(iwph) = wsph(iwph) / 1000.d0
      ENDDO
      wsphmax = wsph(nqstep) 
      CLOSE(iua2ffil)
    ENDIF
    ! first node broadcasts everything to all nodes
    CALL mp_bcast( a2f_iso, ionode_id, inter_pool_comm )
    CALL mp_bcast( wsph, ionode_id, inter_pool_comm )
    CALL mp_bcast( wsphmax, ionode_id, inter_pool_comm )
    CALL mp_barrier(inter_pool_comm)
    !
    WRITE(stdout,'(/5x,a/)') 'Finish reading a2f file '
    !
    IF ( .not. ALLOCATED(memlt_pool) ) ALLOCATE(memlt_pool(npool))
    memlt_pool(:) = 0.d0
    !
    RETURN
    !
    END SUBROUTINE read_a2f
    !
    !-----------------------------------------------------------------------
    SUBROUTINE read_frequencies
    !-----------------------------------------------------------------------
    !
    ! read the frequencies obtained from a previous epw run
    !
    USE io_global, ONLY : stdout, ionode_id
    USE io_epw,    ONLY : iufilfreq
    USE io_files,  ONLY : prefix, tmp_dir
    USE phcom,     ONLY : nmodes
    USE elph2,   ONLY : nqtotf, wf, wqf, xqf
    USE eliashbergcom, ONLY : wsphmax
    USE constants_epw, ONLY : ryd2ev
    USE mp_global, ONLY : inter_pool_comm
    USE mp_world,  ONLY : mpime
    USE mp,        ONLY : mp_bcast, mp_barrier, mp_sum
    !
    IMPLICIT NONE
    !
    INTEGER :: ios
    !! integer variable for I/O control
    INTEGER :: iq
    !! Counter on q points
    INTEGER :: imode
    !! Counter on modes
    CHARACTER (len=256) :: filfreq
    !
    ! read frequencies from file
    IF ( mpime .eq. ionode_id ) THEN
      filfreq = trim(tmp_dir) // trim(prefix) // '.freq'
      !OPEN(iufilfreq, file=filfreq, status='unknown', form='formatted', err=100, iostat=ios)
      OPEN(iufilfreq, file=filfreq, status='unknown', form='unformatted', err=100, iostat=ios)
100   CALL errore('read_frequencies','opening file '//filfreq,abs(ios))
      !READ(iufilfreq,'(2i7)') nqtotf, nmodes
      READ(iufilfreq) nqtotf, nmodes
    ENDIF
    CALL mp_bcast( nqtotf, ionode_id, inter_pool_comm )
    CALL mp_bcast( nmodes, ionode_id, inter_pool_comm )
    !
    IF ( .not. ALLOCATED(wf) )  ALLOCATE(wf(nmodes,nqtotf))
    IF ( .not. ALLOCATED(wqf) ) ALLOCATE(wqf(nqtotf))
    IF ( .not. ALLOCATED(xqf) ) ALLOCATE(xqf(3,nqtotf))
    wf(:,:) = 0.d0
    wqf(:) = 1.d0 / dble(nqtotf)
    xqf(:,:) = 0.d0
    !
    IF ( mpime .eq. ionode_id ) THEN
      DO iq = 1, nqtotf ! loop over q-points
         !READ(iufilfreq,'(3f15.9)') xqf(1,iq), xqf(2,iq), xqf(3,iq)
         !READ(iufilfreq,'(20ES20.10)') (wf(imode,iq), imode=1,nmodes)
         READ(iufilfreq) xqf(1,iq), xqf(2,iq), xqf(3,iq)
         DO imode = 1, nmodes
            READ(iufilfreq) wf(imode,iq)
         ENDDO
      ENDDO 
      CLOSE(iufilfreq)
      ! go from Ryd to eV
      wf(:,:) = wf(:,:) * ryd2ev ! in eV
      wsphmax = 1.1d0 * maxval( wf(:,:) ) ! increase by 10%
    ENDIF
    ! first node broadcasts everything to all nodes
    CALL mp_bcast( wf, ionode_id, inter_pool_comm )
    CALL mp_bcast( xqf, ionode_id, inter_pool_comm )
    CALL mp_bcast( wsphmax, ionode_id, inter_pool_comm )
    CALL mp_barrier(inter_pool_comm)
    !
    WRITE(stdout,'(/5x,a/)') 'Finish reading .freq file '
    !
    RETURN
    !
    END SUBROUTINE read_frequencies
    !
    !-----------------------------------------------------------------------
    SUBROUTINE read_eigenvalues
    !-----------------------------------------------------------------------
    !!
    !! read the eigenvalues obtained from a previous epw run
    !!
    USE kinds,         ONLY : DP
    USE io_global,     ONLY : stdout, ionode_id
    USE io_files,      ONLY : prefix, tmp_dir
    USE pwcom,         ONLY : ef
    USE epwcom,        ONLY : nkf1, nkf2, nkf3, degaussw, fsthick, mp_mesh_k
    USE eliashbergcom, ONLY : nkfs, nbndfs, dosef, ef0, ekfs, wkfs, xkfs, w0g
    USE constants_epw, ONLY : ryd2ev
    USE io_epw,        ONLY : iufilegnv
    USE mp_global, ONLY : inter_pool_comm
    USE mp_world,  ONLY : mpime
    USE mp,        ONLY : mp_bcast, mp_barrier, mp_sum
    ! 
    IMPLICIT NONE
    !
    INTEGER :: ios
    !! integer variable for I/O control
    INTEGER :: ik
    !! Counter on k-points
    INTEGER :: ibnd
    !! Counter on bands
    INTEGER ::  nkftot
    !! Number of k-points
    INTEGER ::  n, nbnd_
    !! Band indexes
    !
    REAL(DP), ALLOCATABLE :: ekf_(:,:)
    !! Temporary eigenvalues on the k point grid
    !
    CHARACTER (len=256) :: filegnv
    REAL(DP), EXTERNAL :: w0gauss
    !
    IF ( mpime .eq. ionode_id ) THEN
      !
      ! SP: Needs to be initialized
      nbnd_ = 0 
      nkfs = 0
      !
      ! read eigenvalues on the irreducible fine k-mesh
      !  
      filegnv = trim(tmp_dir) // trim(prefix) // '.egnv'
      !OPEN(iufilegnv, file=filegnv, status='unknown', form='formatted', err=100, iostat=ios)
      OPEN(iufilegnv, file=filegnv, status='unknown', form='unformatted', err=100, iostat=ios)
100   CALL errore('read_eigenvalues','opening file '//filegnv,abs(ios))
      !
      !READ(iufilegnv,'(5i7)') nkftot, nkf1, nkf2, nkf3, nkfs 
      !READ(iufilegnv,'(i7,5ES20.10)') nbnd_, ef, ef0, dosef, degaussw, fsthick
      READ(iufilegnv) nkftot, nkf1, nkf2, nkf3, nkfs
      READ(iufilegnv) nbnd_, ef, ef0, dosef, degaussw, fsthick
      degaussw = degaussw * ryd2ev
      ef0 = ef0 * ryd2ev
      ef = ef * ryd2ev
      fsthick = fsthick * ryd2ev
      dosef = dosef / ryd2ev
      WRITE(stdout,'(5x,a32,ES20.10)') 'Fermi level (eV) = ', ef0
      WRITE(stdout,'(5x,a32,ES20.10)') 'DOS(states/spin/eV/Unit Cell) = ', dosef
      WRITE(stdout,'(5x,a32,ES20.10)') 'Electron smearing (eV) = ', degaussw
      WRITE(stdout,'(5x,a32,ES20.10)') 'Fermi window (eV) = ', fsthick
      IF ( mp_mesh_k) THEN 
         WRITE(stdout,'(5x,a,i9,a,i9)') 'Nr irreducible k-points within the Fermi shell = ', nkfs, ' out of ', nkftot
      ELSE
         WRITE(stdout,'(5x,a,i9,a,i9)') 'Nr k-points within the Fermi shell = ', nkfs, ' out of ', nkftot
      ENDIF
    ENDIF
    ! first node broadcasts everything to all nodes
    CALL mp_bcast( nkf1, ionode_id, inter_pool_comm )
    CALL mp_bcast( nkf2, ionode_id, inter_pool_comm )
    CALL mp_bcast( nkf3, ionode_id, inter_pool_comm )
    CALL mp_bcast( nkfs, ionode_id, inter_pool_comm )
    CALL mp_bcast( degaussw, ionode_id, inter_pool_comm )
    CALL mp_bcast( ef0, ionode_id, inter_pool_comm )
    CALL mp_bcast( dosef, ionode_id, inter_pool_comm )
    CALL mp_bcast( fsthick, ionode_id, inter_pool_comm )
    CALL mp_bcast( ef, ionode_id, inter_pool_comm )
    !
    IF ( .not. ALLOCATED(wkfs) ) ALLOCATE(wkfs(nkfs))
    IF ( .not. ALLOCATED(xkfs) ) ALLOCATE(xkfs(3,nkfs))
    wkfs(:) = 0.d0
    xkfs(:,:) = 0.d0
    !
    IF ( mpime .eq. ionode_id ) THEN
      !
      ! at each k-point keep only the bands within the Fermi shell
      !
      ALLOCATE(ekf_(nbnd_,nkfs))
      ekf_(:,:) = 0.d0
      !
      ! nbndfs - nr of bands within the Fermi shell
      !
      nbndfs = 0
      DO ik = 1, nkfs ! loop over irreducible k-points
         !READ(iufilegnv,'(4f15.9)') wkfs(ik), xkfs(1,ik), xkfs(2,ik), xkfs(3,ik)
         READ(iufilegnv) wkfs(ik), xkfs(1,ik), xkfs(2,ik), xkfs(3,ik)
         DO ibnd = 1, nbnd_
            !READ(iufilegnv,'(ES20.10)') ekf_(ibnd,ik)
            READ(iufilegnv) ekf_(ibnd,ik)
         ENDDO
         n = 0
         DO ibnd = 1, nbnd_
            ! go from Ryd to eV
            ekf_(ibnd,ik) = ekf_(ibnd,ik) * ryd2ev
            IF ( abs( ekf_(ibnd,ik) - ef0 ) .lt. fsthick ) THEN
               n = n + 1
               IF ( nbndfs .lt. n ) nbndfs = n
            ENDIF
         ENDDO
      ENDDO
      WRITE(stdout,'(5x,i7,a/)') nbndfs, ' bands within the Fermi window'
      CLOSE(iufilegnv)
      ! 
    ENDIF
    ! first node broadcasts everything to all nodes
    CALL mp_bcast( nbndfs, ionode_id, inter_pool_comm )
    CALL mp_bcast( wkfs, ionode_id, inter_pool_comm )
    CALL mp_bcast( xkfs, ionode_id, inter_pool_comm )
    CALL mp_barrier(inter_pool_comm)
    !
    IF ( .not. ALLOCATED(ekfs) ) ALLOCATE(ekfs(nbndfs,nkfs))
    IF ( .not. ALLOCATED(w0g) )  ALLOCATE(w0g(nbndfs,nkfs))
    ! sanity choice
    ekfs(:,:) = ef0 - 10.d0 * fsthick
    w0g(:,:) = 0.d0
    IF ( mpime .eq. ionode_id ) THEN
      DO ik = 1, nkfs ! loop over k-points
         n = 0
         DO ibnd = 1, nbnd_
            IF ( abs( ekf_(ibnd,ik) - ef0 ) .lt. fsthick ) THEN
               n = n + 1
               ekfs(n,ik) = ekf_(ibnd,ik)
               w0g(n,ik) = w0gauss( ( ekfs(n,ik) - ef0 ) / degaussw, 0 ) / degaussw
            ENDIF
         ENDDO
      ENDDO
      IF ( ALLOCATED(ekf_) ) DEALLOCATE(ekf_)
    ENDIF
    ! first node broadcasts everything to all nodes
    CALL mp_bcast( ekfs, ionode_id, inter_pool_comm )
    CALL mp_bcast( w0g, ionode_id, inter_pool_comm )
    CALL mp_barrier(inter_pool_comm)
    !
    WRITE(stdout,'(/5x,a/)') 'Finish reading .egnv file '
    !
    RETURN
    !
    END SUBROUTINE read_eigenvalues
    !
    !-----------------------------------------------------------------------
    SUBROUTINE read_kqmap
    !-----------------------------------------------------------------------
    !
    ! read the map index of k+(sign)q on the k-mesh
    !
    USE kinds,     ONLY : DP
    USE io_global, ONLY : stdout, ionode_id
    USE io_epw,    ONLY : iufilikmap
    USE io_files,  ONLY : prefix, tmp_dir
    USE symm_base, ONLY : t_rev, time_reversal, s, set_sym_bl
    USE phcom,     ONLY : nmodes
    USE epwcom,    ONLY : nkf1, nkf2, nkf3, mp_mesh_k
    USE elph2,     ONLY : nqtotf, xqf
    USE eliashbergcom, ONLY : ixkff, xkff, ixkf, xkfs, nkfs, ixkqf, ixqfs, nbndfs, nqfs, memlt_pool
    USE constants_epw, ONLY : eps5
    USE mp_global,     ONLY : npool
    USE symm_base, ONLY : nrot
    USE mp_global, ONLY : inter_pool_comm
    USE mp_world,  ONLY : mpime
    USE mp,        ONLY : mp_bcast, mp_barrier, mp_sum
    ! 
    IMPLICIT NONE
    !
    INTEGER :: i, j, k, ik, nk, n
    !! Counter on k points
    INTEGER :: iq
    !! Counter on q points
    INTEGER :: nkq
    !! Index of k+sign*q on the fine k-mesh
    INTEGER :: nkftot
    !! Total number of k points
    INTEGER :: nkf_mesh
    !! Nr. of k points read from .ikmap file
    INTEGER :: lower_bnd
    !! Lower bounds index after k or q paral
    INTEGER :: upper_bnd
    !! Upper bounds index after k or q paral
    INTEGER :: nks
    !! Number of non-equivalent k points
    INTEGER :: ns
    !! Counter on rotation operations
    INTEGER :: ios
    !! Integer variable for I/O control
    INTEGER :: imelt
    !! Memory allocated
    INTEGER, ALLOCATABLE :: equiv_(:)
    !! Index of equivalence of k points
    INTEGER, ALLOCATABLE :: index_(:,:)
    !! Index of q-point on the full q-mesh for which k+sign*q is within the Fermi shell
    !
    REAL(kind=DP) :: xk(3)
    !! coordinates of k points
    REAL(kind=DP) :: xq(3)
    !! coordinates of q points
    REAL(kind=DP) :: xkr(3)
    !! coordinates of k points
    REAL(DP) :: xx, yy, zz
    !! Temporary variables
    !
    LOGICAL :: in_the_list
    !! Check if k point is in the list
    CHARACTER (len=256) :: filikmap
    !! Name of the file
    !
    IF ( .not. ALLOCATED(memlt_pool) ) ALLOCATE(memlt_pool(npool))
    memlt_pool(:) = 0.d0
    !
    ! get the size of arrays for frequency and eigenvalue variables allocated in 
    ! read_frequencies and read_eigenvalues
    imelt = ( nmodes + 4 ) * nqtotf + ( 4 + 2 * nbndfs ) * nkfs
    CALL mem_size_eliashberg( imelt )
    !
    nkftot = nkf1 * nkf2 * nkf3
    !
    ! get the size of required memory for ixkff  
    imelt = nkftot
    CALL mem_integer_size_eliashberg( imelt )
    !
    IF ( .not. ALLOCATED(ixkff) ) ALLOCATE(ixkff(nkftot))
    ixkff(:) = 0
    !
    IF ( mpime .eq. ionode_id ) THEN
      !
      filikmap = trim(tmp_dir) // trim(prefix) // '.ikmap'
      !OPEN(iufilikmap, file=filikmap, status='old', form='formatted', err=100, iostat=ios)
      OPEN(iufilikmap, file=filikmap, status='old', form='unformatted', err=100, iostat=ios)
100   CALL errore('read_kqmap','opening file '//filikmap,abs(ios))
      !
      ! nkf_mesh - Total number of k points
      !          - These are irreducible k-points if mp_mesh_k = .true.
      READ(iufilikmap) nkf_mesh
      !
      IF ( .not. ALLOCATED(ixkf) ) ALLOCATE(ixkf(nkf_mesh))
      ixkf(:) = 0
      !
      DO ik = 1, nkf_mesh
         !READ(iufilikmap,'(i9)') ixkf(ik)
         READ(iufilikmap) ixkf(ik)
      ENDDO
      CLOSE(iufilikmap)
      !
      IF ( mp_mesh_k ) CALL set_sym_bl( ) 
      !
      IF ( .not. ALLOCATED(xkff) )  ALLOCATE(xkff(3,nkftot))
      xkff(:,:) = 0.d0
      !
      DO i = 1, nkf1
         DO j = 1, nkf2
            DO k = 1, nkf3
               ik = (i-1)*nkf2*nkf3 + (j-1)*nkf3 + k
               xkff(1,ik) = dble(i-1) / dble(nkf1)
               xkff(2,ik) = dble(j-1) / dble(nkf2)
               xkff(3,ik) = dble(k-1) / dble(nkf3)
            ENDDO
         ENDDO
      ENDDO
      !
      IF ( .not. ALLOCATED(equiv_) )  ALLOCATE(equiv_(nkftot))
      !  equiv_(nk) =nk : k-point nk is not equivalent to any previous k-point
      !  equiv_(nk)!=nk : k-point nk is equivalent to k-point equiv(nk)
      !
      DO nk = 1, nkftot
         equiv_(nk) = nk
      ENDDO
      !
      IF ( mp_mesh_k ) THEN 
         DO nk = 1, nkftot
            !  check if this k-point has already been found equivalent to another
            IF ( equiv_(nk) .eq. nk ) THEN
               !  check if there are equivalent k-point to this in the list
               !  (excepted those previously found to be equivalent to another)
               !  check both k and -k
               DO ns = 1, nrot
                  DO i = 1, 3
                     xkr(i) = s(i,1,ns) * xkff(1,nk) &
                            + s(i,2,ns) * xkff(2,nk) &
                            + s(i,3,ns) * xkff(3,nk)
                     xkr(i) = xkr(i) - nint( xkr(i) )
                  ENDDO
                  IF ( t_rev(ns) .eq. 1 ) xkr = -xkr
                  xx = xkr(1)*nkf1
                  yy = xkr(2)*nkf2
                  zz = xkr(3)*nkf3
                  in_the_list = abs( xx-nint(xx) ) .le. eps5 .AND. &
                                abs( yy-nint(yy) ) .le. eps5 .AND. &
                                abs( zz-nint(zz) ) .le. eps5
                  IF ( in_the_list ) THEN
                     i = mod( nint( xkr(1)*nkf1 + 2*nkf1), nkf1 ) + 1
                     j = mod( nint( xkr(2)*nkf2 + 2*nkf2), nkf2 ) + 1
                     k = mod( nint( xkr(3)*nkf3 + 2*nkf3), nkf3 ) + 1
                     n = (k-1) + (j-1)*nkf3 + (i-1)*nkf2*nkf3 + 1
                     IF ( n .gt. nk .AND. equiv_(n) .eq. n ) THEN
                        equiv_(n) = nk
                     ELSE
                        IF ( equiv_(n) .ne. nk .OR. n .lt. nk ) CALL errore('kmesh_fine', &
                           'something wrong in the checking algorithm',1)
                     ENDIF
                  ENDIF
                  IF ( time_reversal ) THEN
                     xx = -xkr(1)*nkf1
                     yy = -xkr(2)*nkf2
                     zz = -xkr(3)*nkf3
                     in_the_list = abs( xx-nint(xx) ) .le. eps5 .AND. &
                                   abs( yy-nint(yy) ) .le. eps5 .AND. &
                                   abs( zz-nint(zz) ) .le. eps5
                     IF ( in_the_list ) THEN
                        i = mod( nint( -xkr(1)*nkf1 + 2*nkf1), nkf1 ) + 1
                        j = mod( nint( -xkr(2)*nkf2 + 2*nkf2), nkf2 ) + 1
                        k = mod( nint( -xkr(3)*nkf3 + 2*nkf3), nkf3 ) + 1
                        n = (k-1) + (j-1)*nkf3 + (i-1)*nkf2*nkf3 + 1
                        IF ( n .gt. nk .AND. equiv_(n) .eq. n ) THEN
                           equiv_(n) = nk
                        ELSE
                           IF ( equiv_(n) .ne. nk .OR. n .lt. nk ) CALL errore('kmesh_fine', &
                              'something wrong in the checking algorithm',2)
                        ENDIF
                     ENDIF
                  ENDIF
               ENDDO
            ENDIF
         ENDDO
      ENDIF
      !
      !  define index of k on the full mesh (ixkff) using index of k-point within the
      !  Fermi shell (ixkf)
      !
      nks = 0
      DO nk = 1, nkftot
         IF ( equiv_(nk) .eq. nk ) THEN
            nks = nks + 1
            ixkff(nk) = ixkf(nks)
         ELSE
            ixkff(nk) = ixkff(equiv_(nk))
         ENDIF
      ENDDO
      IF ( nks .ne. nkf_mesh) CALL errore('read_kmap_mp', 'something wrong with the mesh',1)
      !
      IF ( ALLOCATED(equiv_) ) DEALLOCATE(equiv_)
      IF ( ALLOCATED(ixkf) )   DEALLOCATE(ixkf)
      !
    ENDIF
    CALL mp_bcast( ixkff, ionode_id, inter_pool_comm )
    CALL mp_barrier(inter_pool_comm)
    !
    CALL fkbounds( nkfs, lower_bnd, upper_bnd )
    !
    ! get the size of required memory for ixkqf, nqfs, index_
    imelt = ( nqtotf + 1 ) * nkfs + ( upper_bnd - lower_bnd + 1 ) * nqtotf
    CALL mem_integer_size_eliashberg( imelt )
    !
    IF ( .not. ALLOCATED(ixkqf) ) ALLOCATE(ixkqf(nkfs,nqtotf))
    IF ( .not. ALLOCATED(nqfs) )  ALLOCATE(nqfs(nkfs))
    IF ( .not. ALLOCATED(index_) ) ALLOCATE(index_(lower_bnd:upper_bnd,nqtotf))
    ixkqf(:,:) = 0
    nqfs(:) = 0
    index_(:,:) = 0
    !
    !
    ! find the index of k+sign*q on the fine k-mesh
    ! nkfs - total nr. of k-points within the Fermi shell (fine mesh)
    !      - these are irreducible k-points if mp_mesh_k=.true.
    ! nqtotf - total nr of q-points on the fine mesh
    !
    DO ik = lower_bnd, upper_bnd
       DO iq = 1, nqtotf
          xk(:) = xkfs(:,ik)
          xq(:) = xqf(:,iq)
          !
          !  nkq - index of k+sign*q on the full fine k-mesh.
          !
          CALL kpmq_map( xk, xq, +1, nkq )
          !
          !  ixkqf(ik,iq) - index of k+sign*q on the fine k-mesh within the Fermi shell
          !
          ixkqf(ik,iq) = ixkff(nkq)
          !
          ! nqfs(ik) - nr of q-points at each k-point for which k+sign*q is within the Fermi shell 
          ! index_   - index q-point on the full q-mesh for which k+sign*q is within the Fermi shell
          !
          IF ( ixkqf(ik,iq) .gt. 0 ) THEN
             nqfs(ik) = nqfs(ik) + 1
             index_(ik,nqfs(ik)) = iq
          ENDIF
       ENDDO
    ENDDO
    !
    ! collect contributions from all pools (sum over k-points)
    CALL mp_sum( ixkqf, inter_pool_comm )
    CALL mp_sum( nqfs,  inter_pool_comm )
    CALL mp_barrier(inter_pool_comm)
    !
    ! get the size of required memory for ixqfs
    imelt = nkfs * maxval(nqfs(:))
    CALL mem_integer_size_eliashberg( imelt )
    !
    IF ( .not. ALLOCATED(ixqfs) ) ALLOCATE(ixqfs(nkfs,maxval(nqfs(:))))
    ixqfs(:,:) = 0
    !
    DO ik = lower_bnd, upper_bnd
       DO iq = 1, nqfs(ik)
          !
          ! ixqfs - index q-point on the full q-mesh for which k+sign*q is within the Fermi shell 
          !
          ixqfs(ik,iq) = index_(ik,iq)
       ENDDO
    ENDDO
    !
    ! collect contributions from all pools (sum over k-points)
    CALL mp_sum( ixqfs, inter_pool_comm )
    CALL mp_barrier(inter_pool_comm)
    !
    IF ( ALLOCATED(index_) ) DEALLOCATE(index_)
    IF ( ALLOCATED(xqf) )    DEALLOCATE(xqf)
    !
    ! remove memory allocated for index_
    imelt = nqtotf * ( upper_bnd - lower_bnd + 1 ) 
    CALL mem_integer_size_eliashberg( -imelt )
    !
    ! remove memory allocated for xqf
    imelt = 3 * nqtotf
    CALL mem_size_eliashberg( -imelt )
    !
    WRITE(stdout,'(/5x,a,i9/)') 'Max nr of q-points = ', maxval(nqfs(:))  
    WRITE(stdout,'(/5x,a/)') 'Finish reading .ikmap files'
    !
    RETURN
    !
    END SUBROUTINE read_kqmap
    !
    !-----------------------------------------------------------------------
    SUBROUTINE read_ephmat
    !-----------------------------------------------------------------------
    !!
    !! Read the electron-phonon matrix elements 
    !!
    USE kinds,         ONLY : DP
    USE io_global,     ONLY : stdout
    USE io_epw,        ONLY : iufileph
    USE io_files,      ONLY : prefix, tmp_dir
    USE phcom,         ONLY : nmodes
    USE elph2,         ONLY : nqtotf, wf
    USE epwcom,        ONLY : eps_acustic, fsthick
    USE eliashbergcom, ONLY : nkfs, nbndfs, ef0, ekfs, g2, ixkqf, nqfs
    USE constants_epw, ONLY : ryd2ev
    USE mp,            ONLY : mp_barrier, mp_bcast, mp_sum
    USE mp_global,     ONLY : inter_pool_comm, npool
    !  
    IMPLICIT NONE
    !
    INTEGER :: ios
    !! integer variable for I/O control
    INTEGER :: ik
    !! Counter on k-points
    INTEGER :: iq
    !! Counter on q-points
    INTEGER :: ibnd
    !! Counter on bands
    INTEGER :: jbnd
    !! Counter on bands
    INTEGER :: imode
    !! Counter on modes
    INTEGER :: nnk
    !! Number of k-points within the Fermi shell
    INTEGER :: nnq(nkfs)
    !! Number of k+q points within the Fermi shell for a given k-point
    INTEGER :: ipool
    !! Counter on pools
    INTEGER :: lower_bnd
    !! Lower bounds index after k or q paral
    INTEGER :: upper_bnd
    !! Upper bounds index after k or q paral
    INTEGER :: tmp_pool_id
    !! Pool index read from file
    INTEGER :: nkpool(npool)
    !! nkpool(ipool) - sum of nr. of k points from pool 1 to pool ipool
    INTEGER :: nmin
    !! Lower bound index for .ephmat file read in current pool
    INTEGER :: nmax
    !! Lower bound index for .ephmat file read in current pool
    INTEGER :: nks
    !! Counter on k points within the Fermi shell
    INTEGER :: imelt
    !! Memory allocated
    !
    REAL(DP) :: gmat
    !! Electron-phonon matrix element square
    !
    CHARACTER (len=256) :: filephmat
    CHARACTER (len=3) :: filelab
    !
    CALL fkbounds( nkfs, lower_bnd, upper_bnd )
    !
    ! get the size of the e-ph matrices that need to be stored in each pool
    imelt = ( upper_bnd - lower_bnd + 1 ) * maxval(nqfs(:)) * nbndfs**2 * nmodes
    CALL mem_size_eliashberg( imelt ) 
    !
    IF ( .not. ALLOCATED(g2) ) ALLOCATE(g2(lower_bnd:upper_bnd,maxval(nqfs(:)),nbndfs,nbndfs,nmodes))
    g2(:,:,:,:,:) = 0.d0
    !
    ! go from Ryd to eV
    ! eps_acustic is given in units of cm-1 in the input file and converted to Ryd in epw_readin
    eps_acustic = eps_acustic * ryd2ev
    !
    WRITE(stdout,'(/5x,a/)') 'Start reading .ephmat files'
    !
    DO ipool = 1, npool ! nr of pools 
       CALL set_ndnmbr(0,ipool,1,npool,filelab)
#if defined(__MPI)
       filephmat = trim(tmp_dir) // trim(prefix) // '.ephmat' // filelab
#else
       filephmat = trim(tmp_dir) // trim(prefix) // '.ephmat'
#endif
       !OPEN(iufileph, file=filephmat, status='old', form='formatted', err=100, iostat=ios)
       OPEN(iufileph, file=filephmat, status='old', form='unformatted', err=100, iostat=ios)
100 CALL errore('read_ephmat','opening file '//filephmat,abs(ios))
       !READ(iufileph,'(2i7)') tmp_pool_id, nkpool(ipool)
       READ(iufileph) tmp_pool_id, nkpool(ipool)
       IF ( ipool .ne. tmp_pool_id )  CALL errore('read_ephmat', &
           'npool should be equal to the number of .ephmat files',1)
       IF ( ipool .gt. 1 ) & 
          nkpool(ipool) = nkpool(ipool) + nkpool(ipool-1)
       !WRITE(stdout,'(2i7)') tmp_pool_id, nkpool(ipool)
       CLOSE(iufileph)
    ENDDO
    CALL mp_barrier(inter_pool_comm)
    !
    ! since the nkfs k-points within the Fermi shell are not evenly distrubed
    ! among the .ephmat files, we re-distribute them here among the npool-pools
    nmin = npool
    nmax = npool
    DO ipool = npool, 1, -1
       IF ( lower_bnd .le. nkpool(ipool) ) THEN
          nmin = ipool
       ENDIF
       IF ( upper_bnd .le. nkpool(ipool) ) THEN
          nmax = ipool
       ENDIF
    ENDDO
    !
    nnk = 0
    nnq(:) = 0
    DO ipool = 1, npool ! nr of pools 
       CALL set_ndnmbr(0,ipool,1,npool,filelab)
#if defined(__MPI)
       filephmat = trim(tmp_dir) // trim(prefix) // '.ephmat' // filelab
#else
       filephmat = trim(tmp_dir) // trim(prefix) // '.ephmat'
#endif     
       OPEN(iufileph, file=filephmat, status='old', form='unformatted')
       READ(iufileph) tmp_pool_id, nks
       IF ( ipool .ge. nmin .AND. ipool .le. nmax ) THEN
          DO iq = 1, nqtotf ! loop over q-points 
             DO ik = 1, nks ! loop over k-points in the pool
                IF ( ixkqf(ik+nnk,iq) .gt. 0 ) THEN 
                   nnq(ik+nnk) = nnq(ik+nnk) + 1
                   DO imode = 1, nmodes ! loop over phonon modes
                      DO ibnd = 1, nbndfs ! loop over iband's 
                         IF ( abs( ekfs(ibnd,ik+nnk) - ef0 ) .lt. fsthick ) THEN
                            DO jbnd = 1, nbndfs ! loop over jband's 
                               IF ( abs( ekfs(jbnd,ixkqf(ik+nnk,iq)) - ef0 ) .lt. fsthick ) THEN
                                  !READ(iufileph,'(ES20.10)') gmat
                                  READ(iufileph) gmat
                                  IF ( ik+nnk .ge. lower_bnd .AND. ik+nnk .le. upper_bnd ) THEN
                                     ! go from Ryd to eV
                                     IF ( wf(imode,iq) .gt. eps_acustic ) THEN
                                        g2(ik+nnk,nnq(ik+nnk),ibnd,jbnd,imode) = gmat * ryd2ev * ryd2ev
                                     ELSE
                                        g2(ik+nnk,nnq(ik+nnk),ibnd,jbnd,imode) = 0.0d0
                                     ENDIF
                                  ENDIF
                               ENDIF ! ekq
                            ENDDO ! jbnd
                         ENDIF ! ekk
                      ENDDO ! ibnd
                   ENDDO ! imode
                ENDIF ! ekk and ekq
             ENDDO ! ik
          ENDDO ! iq
          CLOSE(iufileph)
       ENDIF ! ipool
       nnk = nnk + nks
       IF ( ipool .eq. npool .AND. nnk .ne. nkfs )  CALL errore('read_ephmat', &
           'nnk should be equal to nkfs',1)
    ENDDO ! ipool
    !
    CALL mp_barrier(inter_pool_comm)
    !
    WRITE(stdout,'(/5x,a/)') 'Finish reading .ephmat files '
    !
    RETURN
    !
    END SUBROUTINE read_ephmat
    !
    !-----------------------------------------------------------------------
    !                                                                            
    !-----------------------------------------------------------------------
    SUBROUTINE write_ephmat( iq )
    !-----------------------------------------------------------------------
    !!
    !!  This subroutine writes the elph matrix elements in a format required 
    !!  by Eliashberg equations
    !! 
    !!  Use matrix elements, electronic eigenvalues and phonon frequencies
    !!  from ep-wannier interpolation
    !!
    !-----------------------------------------------------------------------
    USE kinds,      ONLY : DP
    USE io_global,  ONLY : stdout
    USE io_epw,     ONLY : iufilfreq, iufilegnv, iufileph
    USE io_files,   ONLY : prefix, tmp_dir
    USE phcom,      ONLY : nmodes
    USE epwcom,     ONLY : nbndsub, fsthick, ngaussw, degaussw, & 
                           nkf1, nkf2, nkf3, &
                           efermi_read, fermi_energy
    USE pwcom,      ONLY : ef 
    USE elph2,      ONLY : etf, ibndmin, ibndmax, nkqf, epf17, wkf, nkf, &
                           nqtotf, wf, xqf, nkqtotf, efnew 
    USE eliashbergcom, ONLY : equivk, nkfs, ekfs, wkfs, xkfs, dosef, ixkf, ixkqf, nbndfs
    USE constants_epw, ONLY : ryd2ev, two
    USE mp,         ONLY : mp_barrier, mp_sum
    USE mp_global,  ONLY : inter_pool_comm, my_pool_id, npool
    !
    IMPLICIT NONE
    ! 
    INTEGER, INTENT (in) :: iq
    !! Current q-points
    !
    ! Local variables
    !
    INTEGER :: ik
    !! Counter on the k-point index 
    INTEGER :: ikk
    !! k-point index
    INTEGER :: ikq
    !! q-point index 
    INTEGER :: ibnd
    !! Counter on bands
    INTEGER :: jbnd
    !! Counter on bands
    INTEGER :: imode
    !! Counter on mode
    INTEGER :: fermicount
    !! Number of states on the Fermi surface
    INTEGER :: nkftot
    !! Total number of k+q points 
    INTEGER :: lower_bnd
    !! Lower bounds index after k or q paral
    INTEGER :: upper_bnd
    !! Upper bounds index after k or q paral
    INTEGER :: nks
    !! Number of k-point on the current pool
    INTEGER :: imelt
    !! Memory allocated
    !
    REAL(kind=DP) :: ef0
    !! Fermi energy level
    REAL(kind=DP) :: wq
    !! phonon freq
    REAL(kind=DP):: g2
    !! Electron-phonon matrix element square
    REAL(kind=DP), EXTERNAL :: efermig, dos_ef
    !
    CHARACTER (len=256) :: filfreq, filegnv, filephmat
    CHARACTER (len=3) :: filelab
    !
    ! write phonon frequencies to file
    IF ( my_pool_id == 0 ) THEN
      filfreq = trim(tmp_dir) // trim(prefix) // '.freq' 
      IF ( iq .eq. 1 ) THEN
        OPEN(iufilfreq, file = filfreq, form = 'unformatted')
        WRITE(iufilfreq) nqtotf, nmodes
        WRITE(iufilfreq) xqf(1,iq), xqf(2,iq), xqf(3,iq)
        DO imode = 1, nmodes
           WRITE(iufilfreq) wf(imode,iq)
        ENDDO
        CLOSE(iufilfreq)
      ELSE
        OPEN(iufilfreq, file = filfreq, position='append', form = 'unformatted')
        WRITE(iufilfreq) xqf(1,iq), xqf(2,iq), xqf(3,iq)
        DO imode = 1, nmodes
           WRITE(iufilfreq) wf(imode,iq)
        ENDDO
        CLOSE(iufilfreq)
      ENDIF
    ENDIF
    ! 
    ! Fermi level and corresponding DOS
    !  
    ! since wkf(:,ikq) = 0 these bands do not bring any contribution to ef0 or dosef
    ! 
    IF ( efermi_read ) THEN
      ef0 = fermi_energy 
    ELSE
      ef0 = efnew 
      !ef0 = efermig(etf, nbndsub, nkqf, nelec, wkf, degaussw, ngaussw, 0, isk)
      ! if some bands are skipped (nbndskip.neq.0), nelec has already been recalculated 
      ! in ephwann_shuffle
    ENDIF
    !     
    dosef = dos_ef(ngaussw, degaussw, ef0, etf, wkf, nkqf, nbndsub)
    ! N(Ef) in the equation for lambda is the DOS per spin
    dosef = dosef / two
    !
    ! find the bounds of k-dependent arrays in the parallel case
    nkftot = nkqtotf / 2 
    CALL fkbounds( nkftot, lower_bnd, upper_bnd )
    !
    IF (iq.eq.1) THEN
      !
      ! find fermicount - nr of k-points within the Fermi shell per pool
      ! for mp_mesh_k=true. femicount is the nr of irreducible k-points within the Fermi shell per pool
      ! 
      fermicount = 0
      DO ik = 1, nkf
        !
        ikk = 2 * ik - 1
        ikq = ikk + 1
        !
        IF ( equivk(lower_bnd+ik-1) .eq. lower_bnd+ik-1 ) THEN
          IF ( minval( abs( etf(:,ikk) - ef  ) ) .lt. fsthick ) THEN
             fermicount = fermicount + 1
          ENDIF
        ENDIF
        !
      ENDDO
      !
      ! nks = irr nr of k-points within the Fermi shell (fine mesh)
      nks = fermicount
      !
      ! collect contributions from all pools (sum over k-points)
      CALL mp_sum( nks, inter_pool_comm )
      CALL mp_barrier(inter_pool_comm)
      !
      ! write eigenvalues to file
      IF ( my_pool_id == 0 ) THEN
        filegnv = trim(tmp_dir) // trim(prefix) // '.egnv'
        !OPEN(iufilegnv, file = filegnv, form = 'formatted')
        OPEN(iufilegnv, file = filegnv, form = 'unformatted')
        IF ( nks .ne. nkfs ) CALL errore('write_ephmat', &
          'nks should be equal to nr. of irreducible k-points within the Fermi shell on the fine mesh',1)
        !WRITE(iufilegnv,'(5i7)') nkftot, nkf1, nkf2, nkf3, nks
        !WRITE(iufilegnv,'(i7,5ES20.10)') ibndmax-ibndmin+1, ef, ef0, dosef, degaussw, fsthick
        WRITE(iufilegnv) nkftot, nkf1, nkf2, nkf3, nks
        WRITE(iufilegnv) ibndmax-ibndmin+1, ef, ef0, dosef, degaussw, fsthick
        DO ik = 1, nks
           !WRITE(iufilegnv,'(4f15.9)') wkfs(ik), xkfs(1,ik), xkfs(2,ik), xkfs(3,ik) 
           WRITE(iufilegnv) wkfs(ik), xkfs(1,ik), xkfs(2,ik), xkfs(3,ik) 
           DO ibnd = 1, ibndmax-ibndmin+1
              !WRITE(iufilegnv,'(ES20.10)') ekfs(ibnd,ik)
              WRITE(iufilegnv) ekfs(ibnd,ik)
           ENDDO
        ENDDO
        CLOSE(iufilegnv)
      ENDIF
      !
    ENDIF ! iq
    !
    ! write the e-ph matrix elements in the Bloch representation on the fine mesh
    ! in .ephmat files (one for each pool)
    !
#if defined(__MPI)
    CALL set_ndnmbr(0,my_pool_id+1,1,npool,filelab)
    filephmat = trim(tmp_dir) // trim(prefix) // '.ephmat' // filelab
#else
    filephmat = trim(tmp_dir) // trim(prefix) // '.ephmat'
#endif
    IF ( iq .eq. 1 ) THEN 
       OPEN(iufileph, file = filephmat, form = 'unformatted')
    ELSE
       OPEN(iufileph, file = filephmat, position='append', form = 'unformatted')
    ENDIF
    !
    !IF ( iq .eq. 1 ) WRITE(iufileph,'(2i7)') my_pool_id+1, fermicount
    IF ( iq .eq. 1 ) WRITE(iufileph) my_pool_id+1, fermicount
    !
    ! nkf - nr of k-points in the pool (fine mesh)
    ! for mp_mesh_k = true nkf is nr of irreducible k-points in the pool 
    !
    DO ik = 1, nkf
      !  
      ikk = 2 * ik - 1
      ikq = ikk + 1
      !
      ! go only over irreducible k-points
      !
      IF ( equivk(lower_bnd+ik-1) .eq. (lower_bnd+ik-1) ) THEN 
        !
        ! here we must have ef, not ef0, to be consistent with ephwann_shuffle
        !
        !   IF ( ixkf(lower_bnd+ik-1) .gt. 0 .AND. ixkqf(ixkf(lower_bnd+ik-1),iq) .gt. 0 ) THEN
        ! FG: here it can happen that ixkf is 0 and this leads to ixqf(0,iq) after .and.
        !     modified to prevent crash
        IF ( ixkf(lower_bnd+ik-1) > 0 ) THEN
          IF ( ixkqf(ixkf(lower_bnd+ik-1),iq) > 0 ) THEN
            !
            ! 
            DO imode = 1, nmodes ! phonon modes
              wq = wf(imode, iq)
              !
              DO ibnd = 1, ibndmax-ibndmin+1
                IF ( abs( ekfs(ibnd,ixkf(lower_bnd+ik-1)) - ef0 ) < fsthick ) THEN
                  DO jbnd = 1, ibndmax-ibndmin+1
                    IF ( abs( ekfs(jbnd,ixkqf(ixkf(lower_bnd+ik-1),iq)) - ef0 ) < fsthick ) THEN
                      !
                      ! here we take into account the zero-point sqrt(hbar/2M\omega)
                      ! with hbar = 1 and M already contained in the eigenmodes
                      ! g2 is Ry^2, wkf must already account for the spin factor
                      !
                      g2 = abs( epf17(jbnd, ibnd, imode, ik) )**two / ( two * wq )
                      WRITE(iufileph) g2
                    ENDIF
                  ENDDO ! jbnd
                ENDIF
              ENDDO ! ibnd
            ENDDO ! imode
            !
          ENDIF
        ENDIF ! fsthick
        !
      ENDIF ! irr k-points
    ENDDO ! ik's
    CLOSE(iufileph)
    !
    IF ( iq .eq. nqtotf ) THEN 
       IF ( ALLOCATED(ekfs) )   DEALLOCATE(ekfs)
       IF ( ALLOCATED(wkfs) )   DEALLOCATE(wkfs)
       IF ( ALLOCATED(xkfs) )   DEALLOCATE(xkfs)
       IF ( ALLOCATED(ixkqf) )  DEALLOCATE(ixkqf)
       IF ( ALLOCATED(equivk) ) DEALLOCATE(equivk)
       IF ( ALLOCATED(ixkf) )   DEALLOCATE(ixkf)
       !
       ! remove memory allocated for ekfs, wkfs, xkfs 
       imelt = ( nbndfs + 4 ) * nkfs
       CALL mem_size_eliashberg( -imelt )
       !
       ! remove memory allocated for ixkqf 
       imelt = nqtotf * nkfs
       CALL mem_integer_size_eliashberg( -imelt )
       !
       ! remove memory allocated for equivk, ixkf
       imelt = 2 * nkftot
       CALL mem_integer_size_eliashberg( -imelt )
       !
       WRITE(stdout,'(5x,a32,d24.15)') 'Fermi level (eV) = ', ef0 * ryd2ev
       WRITE(stdout,'(5x,a32,d24.15)') 'DOS(states/spin/eV/Unit Cell) = ', dosef / ryd2ev
       WRITE(stdout,'(5x,a32,d24.15)') 'Electron smearing (eV) = ', degaussw * ryd2ev
       WRITE(stdout,'(5x,a32,d24.15)') 'Fermi window (eV) = ', fsthick * ryd2ev
       WRITE(stdout,'(5x,a)')          ' '
       WRITE(stdout,'(5x,a)')          'Finished writing .ephmat files'
       !
    ENDIF
    !
    RETURN
    !
    END SUBROUTINE write_ephmat
    !                                                                            
    !-----------------------------------------------------------------------
    SUBROUTINE count_kpoints( iq )
    !-----------------------------------------------------------------------
    USE kinds,     ONLY : DP
    USE io_global, ONLY : stdout
    USE epwcom,    ONLY : nbndsub, fsthick, ngaussw, degaussw, & 
                          efermi_read, fermi_energy, mp_mesh_k
    USE pwcom,     ONLY : nelec, ef, isk
    USE elph2,     ONLY : etf, nkqf, wkf, nkf, nkqtotf
    USE constants_epw, ONLY : two
    USE mp,        ONLY : mp_barrier, mp_sum
    USE mp_global, ONLY : inter_pool_comm
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT (in) :: iq
    !! Current q-points
    !
    ! Local variables
    !
    INTEGER :: ik
    !! Counter on the k-point index
    INTEGER :: ikk
    !! k-point index
    INTEGER :: ikq
    !! q-point index
    INTEGER :: fermicount
    !! Number of states on the Fermi surface
    INTEGER :: nks
    !! Number of k-point on the current pool
    !
    REAL(kind=DP) :: ef0
    !! Fermi energy level
    REAL(kind=DP) :: dosef
    !! density of states at the Fermi level
    !
    REAL(DP), EXTERNAL :: efermig, dos_ef
    ! 
    IF (iq.eq.1) THEN
       ! 
       ! Fermi level and corresponding DOS
       !  
       ! since wkf(:,ikq) = 0 these bands do not bring any contribution to ef0 or dosef
       !
       IF ( efermi_read ) THEN
          ef0 = fermi_energy 
       ELSE
          ef0 = efermig(etf, nbndsub, nkqf, nelec, wkf, degaussw, ngaussw, 0, isk)
       ENDIF  
       !     
       dosef = dos_ef(ngaussw, degaussw, ef0, etf, wkf, nkqf, nbndsub)
       ! N(Ef) in the equation for lambda is the DOS per spin
       dosef = dosef / two
       !
       ! fermicount = nr of k-points within the Fermi shell per pool
       !
       fermicount = 0
       DO ik = 1, nkf
          !
          ikk = 2 * ik - 1
          ikq = ikk + 1
          !
          IF ( minval( abs( etf(:,ikk) - ef  ) ) .lt. fsthick ) &
             fermicount = fermicount + 1 
          !
       ENDDO
       !
       ! nks =  nr of k-points within the Fermi shell (fine mesh)
       nks = fermicount
       !
       ! collect contributions from all pools (sum over k-points)
       CALL mp_sum( nks, inter_pool_comm )
       CALL mp_barrier(inter_pool_comm)
       !
       IF ( mp_mesh_k) THEN
          WRITE(stdout,'(5x,a,i9,a,i9)') 'Nr irreducible k-points within the Fermi shell = ', nks, ' out of ', nkqtotf / 2
       ELSE
          WRITE(stdout,'(5x,a,i9,a,i9)') 'Nr k-points within the Fermi shell = ', nks, ' out of ', nkqtotf / 2
       ENDIF
    ENDIF ! iq
    !
    RETURN
    !
    END SUBROUTINE count_kpoints
    !                                                                            
    !-----------------------------------------------------------------------
    !                                                                            
    !-----------------------------------------------------------------------
    SUBROUTINE kmesh_fine
    !-----------------------------------------------------------------------
    !!
    !!   This routine defines the nr. of k-points on the fine k-mesh 
    !!   within the Fermi shell
    !!
    USE kinds,     ONLY : DP
    USE io_global, ONLY : ionode_id, stdout
    USE io_files,  ONLY : prefix, tmp_dir
    USE epwcom,    ONLY : nkf1, nkf2, nkf3, fsthick, mp_mesh_k
    USE pwcom,     ONLY : ef
    USE io_epw,    ONLY : iufilikmap
    USE elph2,     ONLY : xkf, wkf, etf, nkf, nkqtotf, ibndmin, ibndmax
    USE eliashbergcom, ONLY : nkfs, ixkf, equivk, xkfs, wkfs, ekfs, nbndfs, memlt_pool
    USE mp_global, ONLY : inter_pool_comm, npool
    USE mp,        ONLY : mp_bcast, mp_barrier, mp_sum
    USE mp_world,  ONLY : mpime
    !
    IMPLICIT NONE
    !
    INTEGER :: nk
    !! Counter on k points
    INTEGER :: nks
    !! Counter on k points within the Fermi shell
    INTEGER :: ikk
    !! k-point index
    INTEGER :: nkf_mesh
    !! Total number of k points
    !! These are irreducible k-points if mp_mesh_k = .true.
    INTEGER :: lower_bnd
    !! Lower bounds index after k or q paral
    INTEGER :: upper_bnd
    !! Upper bounds index after k or q paral
    INTEGER :: imelt
    !! Memory allocated
    REAL(DP) :: xx, yy, zz
    !!
    REAL(DP), ALLOCATABLE :: xkf_(:,:)
    !! Temporary k point grid
    REAL(DP), ALLOCATABLE :: wkf_(:)
    !! Temporary weights on the k point grid
    REAL(DP), ALLOCATABLE :: ekf_(:,:)
    !! Temporary eigenvalues on the k point grid
    CHARACTER (len=256) :: filikmap
    !! Name of the file
    !
    nkf_mesh = nkqtotf / 2 
    nbndfs = ibndmax - ibndmin + 1
    !
#if defined(__MPI)
    IF ( .not. ALLOCATED(memlt_pool) ) ALLOCATE(memlt_pool(npool))
    memlt_pool(:) = 0.d0
#else
    IF ( .not. ALLOCATED(memlt_pool) ) ALLOCATE(memlt_pool(1))
    memlt_pool(1) = 0.d0
#endif
    !
    ! get the size of required memory for ekf_, wkf_, xkf_
    imelt = ( nbndfs + 4 ) * nkf_mesh 
    CALL mem_size_eliashberg( imelt )
    !
    ! get the size of required memory for ixkf and equivk
    imelt = 2 * nkf_mesh
    CALL mem_integer_size_eliashberg( imelt )
    !
    IF ( .not. ALLOCATED(ekf_) )   ALLOCATE(ekf_(nbndfs,nkf_mesh))
    IF ( .not. ALLOCATED(wkf_) )   ALLOCATE(wkf_(nkf_mesh))
    IF ( .not. ALLOCATED(xkf_) )   ALLOCATE(xkf_(3,nkf_mesh))
    IF ( .not. ALLOCATED(equivk) ) ALLOCATE(equivk(nkf_mesh))
    IF ( .not. ALLOCATED(ixkf) )   ALLOCATE(ixkf(nkf_mesh))
    xkf_(:,:) = 0.d0
    ekf_(:,:) = 0.d0
    wkf_(:) = 0.d0
    equivk(:) = 0
    ixkf(:) = 0
    !
    CALL fkbounds( nkf_mesh, lower_bnd, upper_bnd )
    !
    ! nkf - nr of k-blocks in the pool (fine grid)
    !
    DO nk = 1, nkf
       ikk = 2 * nk - 1
       xkf_(:,lower_bnd+nk-1) = xkf(:,ikk)
       wkf_(lower_bnd+nk-1)   = wkf(ikk)
       ekf_(:,lower_bnd+nk-1) = etf(ibndmin:ibndmax,ikk)
    ENDDO
       !
       ! collect contributions from all pools (sum over k-points)
       CALL mp_sum( ekf_, inter_pool_comm )
       CALL mp_sum( xkf_, inter_pool_comm )
       CALL mp_sum( wkf_, inter_pool_comm )
       CALL mp_barrier(inter_pool_comm)
    !
    IF ( mpime .eq. ionode_id ) THEN
      DO nk = 1, nkf_mesh
         equivk(nk) = nk
      ENDDO
      !
      IF ( mp_mesh_k) THEN
         WRITE(stdout,'(/5x,a,i9/)') 'Nr. of irreducible k-points on the uniform grid: ', nkf_mesh
      ELSE
         WRITE(stdout,'(/5x,a,i9/)') 'Nr. of k-points on the uniform grid: ', nkf_mesh
      ENDIF
      !
      filikmap = trim(tmp_dir) // trim(prefix) // '.ikmap'
      !OPEN(iufilikmap, file = filikmap, form = 'formatted')
      !WRITE(iufilikmap,'(i9)') nkf_mesh
      OPEN(iufilikmap, file = filikmap, form = 'unformatted')
      WRITE(iufilikmap) nkf_mesh
      !
      ! nkfs - find nr of k-points within the Fermi shell (fine grid)
      ! only a fraction of nkf_mesh are contained in the Fermi shell
      !
      ! ixkf - find the index of k-point within the Fermi shell (fine grid)
      ! if the k-point lies outside the Fermi shell the index is 0
      !
      nkfs = 0  
      DO nk = 1, nkf_mesh
         IF ( minval( abs( ekf_(:,nk) - ef  ) ) .lt. fsthick ) THEN
            nkfs = nkfs + 1
            ixkf(nk) = nkfs
         ELSE
            ixkf(nk) = 0
         ENDIF
         !  bring back into to the first BZ
         xx = xkf_(1,nk) * nkf1
         yy = xkf_(2,nk) * nkf2
         zz = xkf_(3,nk) * nkf3
         CALL backtoBZ( xx, yy, zz, nkf1, nkf2, nkf3 )
         xkf_(1,nk) = xx / dble(nkf1)
         xkf_(2,nk) = yy / dble(nkf2)
         xkf_(3,nk) = zz / dble(nkf3)
         !WRITE(iufilikmap,'(i9)') ixkf(nk)
         WRITE(iufilikmap) ixkf(nk)
      ENDDO
      CLOSE(iufilikmap)
      !
    ENDIF
    CALL mp_bcast( nkfs, ionode_id, inter_pool_comm )
    !
    ! get the size of required memory for ekfs, wkfs, xkfs 
    imelt = ( nbndfs + 4 ) * nkfs
    CALL mem_size_eliashberg( imelt )
    ! 
    IF ( .not. ALLOCATED(ekfs) ) ALLOCATE(ekfs(nbndfs,nkfs))
    IF ( .not. ALLOCATED(wkfs) ) ALLOCATE(wkfs(nkfs))
    IF ( .not. ALLOCATED(xkfs) ) ALLOCATE(xkfs(3,nkfs))
    xkfs(:,:) = 0.d0
    wkfs(:) = 0.d0
    ekfs(:,:) = 0.d0
    !
    IF ( mpime .eq. ionode_id ) THEN
      nks = 0
      DO nk = 1, nkf_mesh
         IF ( minval( abs( ekf_(:,nk) - ef  ) ) .lt. fsthick ) THEN
            nks = nks + 1
            IF ( nks .gt. nkf_mesh ) CALL errore('kmesh_fine','too many k-points',1)
            wkfs(nks)   = wkf_(nk)
            xkfs(:,nks) = xkf_(:,nk)
            ekfs(:,nks) = ekf_(:,nk)
         ENDIF
      ENDDO
    ENDIF
    !
    ! first node broadcasts everything to all nodes
    CALL mp_bcast( ixkf, ionode_id, inter_pool_comm )
    CALL mp_bcast( equivk, ionode_id, inter_pool_comm )
    CALL mp_bcast( xkfs, ionode_id, inter_pool_comm )
    CALL mp_bcast( wkfs, ionode_id, inter_pool_comm )
    CALL mp_bcast( ekfs, ionode_id, inter_pool_comm )
    CALL mp_barrier(inter_pool_comm)
    !
    IF ( ALLOCATED(ekf_) ) DEALLOCATE(ekf_)
    IF ( ALLOCATED(xkf_) ) DEALLOCATE(xkf_)
    IF ( ALLOCATED(wkf_) ) DEALLOCATE(wkf_)
    !
    ! remove memory allocated for ekf_, xkf_, wkf_
    imelt = ( nbndfs + 4 ) * nkf_mesh
    CALL mem_size_eliashberg( -imelt )
    !
    WRITE(stdout,'(/5x,a/)') 'Finished writing .ikmap file '
    !
    RETURN
    !
    END SUBROUTINE kmesh_fine
    !
    !-----------------------------------------------------------------------
    SUBROUTINE kqmap_fine
    !-----------------------------------------------------------------------
    !!
    !! this routine finds the index of k+sign*q on the fine k-mesh
    !!
    USE kinds,     ONLY : DP
    USE io_global, ONLY : stdout
    USE symm_base, ONLY : s, t_rev, time_reversal, set_sym_bl
    USE epwcom,    ONLY : nkf1, nkf2, nkf3, mp_mesh_k
    USE elph2,     ONLY : nqtotf, xqf
    USE eliashbergcom, ONLY : ixkff, xkff, ixkf, xkfs, nkfs, ixkqf, ixqfs, nqfs
    USE constants_epw, ONLY : eps5
    USE symm_base, ONLY : nrot
    USE io_global, ONLY : ionode_id
    USE mp_global, ONLY : inter_pool_comm
    USE mp,        ONLY : mp_bcast, mp_barrier, mp_sum
    USE mp_world,  ONLY : mpime
    ! 
    IMPLICIT NONE
    !
    INTEGER :: i, j, k, ik, nk, n
    !! Counter on k points
    INTEGER :: iq
    !! Counter on q points
    INTEGER :: nkq
    !! Index of k+sign*q on the fine k-mesh
    INTEGER :: nkftot
    !! Total number of k points
    INTEGER :: lower_bnd
    !! Lower bounds index after k or q paral
    INTEGER :: upper_bnd
    !! Upper bounds index after k or q paral
    INTEGER :: nks
    !! Number of non-equivalent k points
    INTEGER :: ns
    !! Counter on rotation operations
    INTEGER :: imelt
    !! Memory allocated
    INTEGER, ALLOCATABLE :: equiv_(:)
    !! Index of equivalence of k points
    INTEGER, ALLOCATABLE :: index_(:,:)
    !! Index of q-point on the full q-mesh for which k+sign*q is within the Fermi shell
    REAL(kind=DP) :: xk(3)
    !! coordinates of k points
    REAL(kind=DP) :: xq(3)
    !! coordinates of q points
    REAL(kind=DP) :: xkr(3)
    !! coordinates of k points
    REAL(DP) :: xx, yy, zz
    !! Temporary variables
    LOGICAL :: in_the_list
    !! Check if k point is in the list
    !
    nkftot = nkf1 * nkf2 * nkf3
    !
    ! get the size of required memory for xkff
    imelt = 3 * nkftot 
    CALL mem_size_eliashberg( imelt )
    !
    ! get the size of required memory for ixkff and equiv_
    imelt = 2 * nkftot
    CALL mem_integer_size_eliashberg( imelt )
    !
    IF ( .not. ALLOCATED(xkff) )  ALLOCATE(xkff(3,nkftot))
    IF ( .not. ALLOCATED(ixkff) ) ALLOCATE(ixkff(nkftot))
    xkff(:,:) = 0.d0
    ixkff(:) = 0
    !
    ! to map k+q onto k we need to define the index of k on the full mesh (ixkff) 
    ! using index of the k-point within the Fermi shell (ixkf)
    !
    IF ( mpime .eq. ionode_id ) THEN
      !
      IF ( mp_mesh_k ) CALL set_sym_bl( ) 
      !
      DO i = 1, nkf1
         DO j = 1, nkf2
            DO k = 1, nkf3
               ik = (i-1)*nkf2*nkf3 + (j-1)*nkf3 + k
               xkff(1,ik) = dble(i-1) / dble(nkf1)
               xkff(2,ik) = dble(j-1) / dble(nkf2)
               xkff(3,ik) = dble(k-1) / dble(nkf3)
            ENDDO
         ENDDO
      ENDDO
      !
      IF ( .not. ALLOCATED(equiv_) ) ALLOCATE(equiv_(nkftot))
      !  equiv_(nk) =nk : k-point nk is not equivalent to any previous k-point
      !  equiv_(nk)!=nk : k-point nk is equivalent to k-point equiv(nk)
      !
      DO nk = 1, nkftot
         equiv_(nk) = nk
      ENDDO
      !
      IF ( mp_mesh_k ) THEN
        DO nk = 1, nkftot
          !  check if this k-point has already been found equivalent to another
          IF ( equiv_(nk) .eq. nk ) THEN
            !  check if there are equivalent k-point to this in the list
            !  (excepted those previously found to be equivalent to another)
            !  check both k and -k
            DO ns = 1, nrot
              DO i = 1, 3
                xkr(i) = s(i,1,ns) * xkff(1,nk) &
                       + s(i,2,ns) * xkff(2,nk) &
                       + s(i,3,ns) * xkff(3,nk)
                xkr(i) = xkr(i) - nint( xkr(i) )
              ENDDO
              IF ( t_rev(ns) .eq. 1 ) xkr = -xkr
              xx = xkr(1)*nkf1
              yy = xkr(2)*nkf2
              zz = xkr(3)*nkf3
              in_the_list = abs( xx-nint(xx) ) .le. eps5 .AND. &
                            abs( yy-nint(yy) ) .le. eps5 .AND. &
                            abs( zz-nint(zz) ) .le. eps5
              IF ( in_the_list ) THEN
                i = mod( nint( xkr(1)*nkf1 + 2*nkf1), nkf1 ) + 1
                j = mod( nint( xkr(2)*nkf2 + 2*nkf2), nkf2 ) + 1
                k = mod( nint( xkr(3)*nkf3 + 2*nkf3), nkf3 ) + 1
                n = (k-1) + (j-1)*nkf3 + (i-1)*nkf2*nkf3 + 1
                IF ( n .gt. nk .AND. equiv_(n) .eq. n ) THEN
                   equiv_(n) = nk
                ELSE
                   IF ( equiv_(n) .ne. nk .OR. n .lt. nk ) CALL errore('kmesh_fine', &
                      'something wrong in the checking algorithm',1)
                ENDIF
              ENDIF
              IF ( time_reversal ) THEN
                xx = -xkr(1)*nkf1
                yy = -xkr(2)*nkf2
                zz = -xkr(3)*nkf3
                in_the_list = abs( xx-nint(xx) ) .le. eps5 .AND. &
                              abs( yy-nint(yy) ) .le. eps5 .AND. &
                              abs( zz-nint(zz) ) .le. eps5
                IF ( in_the_list ) THEN
                  i = mod( nint( -xkr(1)*nkf1 + 2*nkf1), nkf1 ) + 1
                  j = mod( nint( -xkr(2)*nkf2 + 2*nkf2), nkf2 ) + 1
                  k = mod( nint( -xkr(3)*nkf3 + 2*nkf3), nkf3 ) + 1
                  n = (k-1) + (j-1)*nkf3 + (i-1)*nkf2*nkf3 + 1
                  IF ( n .gt. nk .AND. equiv_(n) .eq. n ) THEN
                    equiv_(n) = nk
                  ELSE
                    IF ( equiv_(n) .ne. nk .OR. n .lt. nk ) CALL errore('kmesh_fine', &
                       'something wrong in the checking algorithm',2)
                  ENDIF
                ENDIF
              ENDIF
            ENDDO
          ENDIF
        ENDDO
      ENDIF
      !
      ! find index of k on the full mesh (ixkff) using index of k within the Fermi shell (ixkf)
      ! 
      nks = 0
      DO nk = 1, nkftot
        IF ( equiv_(nk) .eq. nk ) THEN
          nks = nks + 1
          ixkff(nk) = ixkf(nks)
        ELSE
          ixkff(nk) = ixkff(equiv_(nk))
        ENDIF
      ENDDO
      !
      IF ( ALLOCATED(equiv_) ) DEALLOCATE(equiv_)
      !
    ENDIF
    CALL mp_bcast( xkff, ionode_id, inter_pool_comm )
    CALL mp_bcast( ixkff, ionode_id, inter_pool_comm )
    CALL mp_barrier(inter_pool_comm)
    !
    IF ( ALLOCATED(xkff) ) DEALLOCATE(xkff)
    !
    ! remove memory allocated for xkff
    imelt = 3 * nkftot
    CALL mem_size_eliashberg( -imelt )
    !
    ! remove memory allocated for equiv_
    imelt = nkftot
    CALL mem_integer_size_eliashberg( -imelt )
    !
    CALL fkbounds( nkfs, lower_bnd, upper_bnd )
    !
    ! get the size of required memory for ixkqf, nqfs, index_
    imelt = ( nqtotf + 1 ) * nkfs + ( upper_bnd - lower_bnd + 1 ) * nqtotf
    CALL mem_integer_size_eliashberg( imelt )
    !
    IF ( .not. ALLOCATED(ixkqf) )  ALLOCATE(ixkqf(nkfs,nqtotf))
    IF ( .not. ALLOCATED(nqfs) )   ALLOCATE(nqfs(nkfs))
    IF ( .not. ALLOCATED(index_) ) ALLOCATE(index_(lower_bnd:upper_bnd,nqtotf))
    ixkqf(:,:) = 0
    nqfs(:) = 0
    index_(:,:) = 0
    !
    ! find the index of k+sign*q on the fine k-mesh
    ! nkfs - total nr. of k-points within the Fermi shell
    !      - these are irreducible k-points if mp_mesh_k = true
    ! nqtotf - total nr of q-points on the fine mesh
    !
    DO ik = lower_bnd, upper_bnd
      DO iq = 1, nqtotf
        xk(:) = xkfs(:,ik)
        xq(:) = xqf(:,iq)
        !
        ! find nkq - index of k+sign*q on the full fine k-mesh.
        !
        CALL kpmq_map( xk, xq, +1, nkq )
        !
        ! find ixkqf(ik,iq) - index of k+sign*q on the fine k-mesh within the Fermi shell
        !
        ixkqf(ik,iq) = ixkff(nkq) 
        !
        ! nqfs(ik) - nr of q-points at each k-point for which k+sign*q is within the Fermi shell 
        ! index_   - index of q-point on the full q-mesh for which k+sign*q is within the Fermi shell
        !
        IF ( ixkqf(ik,iq) .gt. 0 ) THEN
          nqfs(ik) = nqfs(ik) + 1
          index_(ik,nqfs(ik)) = iq
        ENDIF
      ENDDO ! loop over full set of q-points (fine mesh)
    ENDDO ! loop over k-points within the Fermi shell in each pool (fine mesh) 
    !
    ! collect contributions from all pools (sum over k-points)
    CALL mp_sum( ixkqf, inter_pool_comm )
    CALL mp_sum( nqfs,  inter_pool_comm )
    CALL mp_barrier(inter_pool_comm)
    !
    ! get the size of required memory for ixqfs
    imelt = nkfs * maxval(nqfs(:))
    CALL mem_integer_size_eliashberg( imelt )
    !
    IF ( .not. ALLOCATED(ixqfs) ) ALLOCATE(ixqfs(nkfs,maxval(nqfs(:))))
    ixqfs(:,:) = 0
    !
    DO ik = lower_bnd, upper_bnd
      DO iq = 1, nqfs(ik)
        !
        ! ixqfs - index of q-point on the full q-mesh for which k+sign*q is within the Fermi shell 
        !
        ixqfs(ik,iq) = index_(ik,iq)   
      ENDDO
    ENDDO
    !
    ! collect contributions from all pools (sum over k-points)
    CALL mp_sum( ixqfs, inter_pool_comm )
    CALL mp_barrier(inter_pool_comm)
    !
    ! remove memory allocated for ixkff, ixqfs, index_, nqfs
    imelt = nkftot + nkfs * maxval(nqfs(:)) + nqtotf * ( upper_bnd - lower_bnd + 1 ) + nkfs
    CALL mem_integer_size_eliashberg( -imelt )
    !
    IF ( ALLOCATED(ixkff) )  DEALLOCATE(ixkff)
    IF ( ALLOCATED(ixqfs) )  DEALLOCATE(ixqfs)
    IF ( ALLOCATED(index_) ) DEALLOCATE(index_)
    IF ( ALLOCATED(nqfs) )   DEALLOCATE(nqfs)
    !
    IF ( mp_mesh_k) THEN 
      WRITE(stdout,'(/5x,a/)') 'Finished mapping k+sign*q onto the fine irreducibe k-mesh'
    ELSE
      WRITE(stdout,'(/5x,a/)') 'Finished mapping k+sign*q onto the fine k-mesh'
    ENDIF
    ! 
    RETURN
    !
    END SUBROUTINE kqmap_fine
    !
    !-----------------------------------------------------------------------
    SUBROUTINE kpmq_map( xk, xq, sign1, nkq )
    !-----------------------------------------------------------------------
    !!
    !! this routine finds the index of k+q or k-q point on the fine k-mesh
    !!
    USE kinds,     ONLY : DP
    USE epwcom,    ONLY : nkf1, nkf2, nkf3
    USE constants_epw, ONLY : eps5
    USE mp,        ONLY : mp_bcast, mp_barrier
    ! 
    IMPLICIT NONE
    !
    INTEGER, INTENT (in) :: sign1
    !! +1 for searching k+q, -1 for k-q
    INTEGER, INTENT (out) :: nkq
    !! the index of k+sign*q
    ! 
    REAL(kind=DP), INTENT (in) :: xk(3)
    !! coordinates of k points
    REAL(kind=DP), INTENT (in) :: xq(3)
    !! coordinates of q points
    ! 
    ! Local variables
    REAL(DP) :: xx, yy, zz, xxk(3)
    LOGICAL :: in_the_list
    !
    !
    xxk(:) = xk(:) + dble(sign1) * xq(:)
    xx = xxk(1) * nkf1
    yy = xxk(2) * nkf2
    zz = xxk(3) * nkf3
    in_the_list = abs(xx-nint(xx)) .le. eps5 .AND. &
                  abs(yy-nint(yy)) .le. eps5 .AND. &
                  abs(zz-nint(zz)) .le. eps5
    IF ( .not. in_the_list ) CALL errore('kpmq_map','k+q does not fall on k-grid',1)
    !
    !  find the index of this k+q or k-q in the k-grid
    !  make sure xx, yy, zz are in the 1st BZ
    !
    CALL backtoBZ( xx, yy, zz, nkf1, nkf2, nkf3 )
    !
    ! since k- and q- meshes are commensurate, nkq can be easily found
    !
    nkq = nint(xx) * nkf2 * nkf3 + nint(yy) * nkf3 + nint(zz) + 1
    !
    !  Now nkq represents the index of k+sign*q on the fine k-grid.
    !
    RETURN
    ! 
    END SUBROUTINE kpmq_map
    !
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
