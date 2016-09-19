  !                                                                            
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino 
  ! Copyright (C) 2007-2009  Roxana Margine, Feliciano Giustino  
  !                                                                        
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  !----------------------------------------------------------------------- 
  SUBROUTINE eliashberg_init
  !-----------------------------------------------------------------------
  !
  ! This routine initializes the control variables needed to solve the eliashberg 
  ! equations
  !
  USE kinds,         ONLY : DP
  USE io_global,     ONLY : stdout
  USE epwcom,        ONLY : eliashberg, nkf1, nkf2, nkf3, nsiter, nqstep, &
                            nqf1, nqf2, nqf3, ntempxx, nswi, nswfc, nswc, &
                            nstemp, muc, lreal, lpade, liso, limag, laniso, &
                            lacon, kerwrite, kerread, imag_read, fila2f, temps, &
                            wsfc, wscut, tempsmin, tempsmax, rand_q, rand_k
  USE constants_epw, ONLY : pi, kelvin2eV
  USE eliashbergcom, ONLY : estemp, nsw, nsiw, dwsph, wsphmax, wsph
  !
  IMPLICIT NONE
  !
  INTEGER :: itemp, iwph, imelt
  REAL(DP) :: dtemp
  REAL(DP) :: eps=1.0d-6
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
     estemp(:) = temps(:) * kelvin2eV
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
     IF ( ABS(wsfc) < eps .OR. ABS(wscut) < eps .OR. nswfc .eq. 0 .OR. nswc .eq. 0 ) THEN 
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
     IF ( ABS(wscut) < eps ) THEN 
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
    WRITE(iufillambda,'(a12,a30)') '# lambda_nk','  \rho(lambda_nk) scaled to 1.'
    DO ibin = 1, nbink
      WRITE(iufillambda,'(2f21.7)') dbink*dble(ibin), lambda_k_bin(ibin)/maxval(lambda_k_bin(:))
    ENDDO
    CLOSE(iufillambda)
    !
    ! SP: Produced if user really wants it 
    IF ( iverbosity == 2 ) THEN  
      OPEN( unit = iufillambda, file = TRIM(prefix)//".lambda_pairs", form = 'formatted')
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
          WRITE(name1,'(a,a8,i1,a5)') TRIM(prefix),'.lambda_', ibnd, '.cube'
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

