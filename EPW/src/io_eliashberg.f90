  !
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino  
  ! Copyright (C) 2007-2009 Roxana Margine
  ! 
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .
  !
  !----------------------------------------------------------------------
  MODULE io_eliashberg
  !----------------------------------------------------------------------
  !! 
  !! This module contains all the IO part of the superconductivity part of EPW
  !!  
  IMPLICIT NONE
  ! 
  CONTAINS
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
    USE constants_epw, ONLY : kelvin2eV, eps6, zero
    USE io_global, ONLY : ionode_id
    USE mp_global, ONLY : inter_pool_comm
    USE mp_world,  ONLY : mpime
    USE mp,        ONLY : mp_bcast, mp_barrier, mp_sum
    USE superconductivity, ONLY : mem_size_eliashberg, free_energy
    ! 
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: itemp
    !! Counter on temperature
    !
    ! Local variables
    INTEGER :: iw
    !! Counter on frequency
    INTEGER :: ik
    !! Counter on k-point
    INTEGER :: ibnd
    !! Counter on band
    INTEGER :: imelt
    !! Required allocation of memory
    INTEGER :: ios
    !! Status variables when reading a file
    !
    REAL(DP) :: temp
    !! Temperature in K
    REAL(DP) :: eband
    !! Temporary variable for eigenvalue
    REAL(DP) :: omega
    !! Temporary variable for frequency
    REAL(DP) :: weight
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
    gap(:) = zero
    Agap(:,:,:) = zero
    Deltai(:) = zero
    Znormi(:) = zero
    NZnormi(:) = zero
    ADeltai(:,:,:) = zero
    AZnormi(:,:,:) = zero
    NAZnormi(:,:,:) = zero
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
    INTEGER :: iw
    !! Counter on frequency imag-axis
    INTEGER :: ik
    !! Counter on k-points
    INTEGER :: ibnd
    !! Counter on bands
    !
    REAL(DP) :: temp
    !! Temperature in K
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
    SUBROUTINE eliashberg_write_raxis( itemp, cname )
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
    INTEGER, INTENT (in) :: itemp
    !! Counter for temperature
    CHARACTER(len=256), INTENT (in) :: cname
    !
    ! Local variables
    INTEGER :: iw
    !! Counter on frequency real-axis
    INTEGER :: ik
    !! Counter on k-points
    INTEGER :: ibnd
    !! Counter on bands
    !
    REAL(DP) :: temp
    !! Temperature in K
    CHARACTER (len=256) :: name1
    LOGICAL :: lgap
    !! True if gap found
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
    END SUBROUTINE eliashberg_write_raxis
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
    !-----------------------------------------------------------------------
    SUBROUTINE read_a2f
    !-----------------------------------------------------------------------
    !!
    !! Read the eliashberg spectral function from fila2f
    !!
    USE io_global,     ONLY : stdout
    USE epwcom,        ONLY : nqstep, fila2f
    USE eliashbergcom, ONLY : wsphmax, wsph, a2f_iso, memlt_pool
    USE constants_epw, ONLY : zero
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
    a2f_iso(:) = zero
    wsph(:) = zero
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
    USE constants_epw, ONLY : ryd2ev, zero
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
    wf(:,:) = zero
    wqf(:) = 1.d0 / dble(nqtotf)
    xqf(:,:) = zero
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
    USE constants_epw, ONLY : ryd2ev, zero
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
    wkfs(:) = zero
    xkfs(:,:) = zero
    !
    IF ( mpime .eq. ionode_id ) THEN
      !
      ! at each k-point keep only the bands within the Fermi shell
      !
      ALLOCATE(ekf_(nbnd_,nkfs))
      ekf_(:,:) = zero
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
    w0g(:,:) = zero
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
    USE superconductivity, ONLY : mem_size_eliashberg, mem_integer_size_eliashberg
    USE constants_epw, ONLY : eps5, zero
    USE symm_base, ONLY : nrot
    USE mp_global, ONLY : inter_pool_comm, npool
    USE mp_world,  ONLY : mpime
    USE mp,        ONLY : mp_bcast, mp_barrier, mp_sum
    USE division,  ONLY : fkbounds
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
    !! Lower bounds index after k paral
    INTEGER :: upper_bnd
    !! Upper bounds index after k paral
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
    memlt_pool(:) = zero
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
      IF ( .not. ALLOCATED(xkff) )  ALLOCATE(xkff(3,nkftot))
      xkff(:,:) = zero
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
         CALL set_sym_bl( ) 
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
    USE superconductivity, ONLY : mem_size_eliashberg
    USE constants_epw, ONLY : ryd2ev, zero
    USE mp,            ONLY : mp_barrier, mp_bcast, mp_sum
    USE mp_global,     ONLY : inter_pool_comm, npool
    USE division,      ONLY : fkbounds
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
    !! Lower bounds index after k paral
    INTEGER :: upper_bnd
    !! Upper bounds index after k paral
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
    g2(:,:,:,:,:) = zero
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
    SUBROUTINE write_ephmat( iqq, iq, totq )
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
                           wf, xqf, nkqtotf, efnew 
    USE eliashbergcom, ONLY : equivk, nkfs, ekfs, wkfs, xkfs, dosef, ixkf, ixkqf, nbndfs
    USE superconductivity, ONLY : mem_size_eliashberg, mem_integer_size_eliashberg
    USE constants_epw, ONLY : ryd2ev, two
    USE mp,         ONLY : mp_barrier, mp_sum
    USE mp_global,  ONLY : inter_pool_comm, my_pool_id, npool
    USE division,   ONLY : fkbounds
    !
    IMPLICIT NONE
    ! 
    INTEGER, INTENT (in) :: iqq
    !! Current q-points from selecq
    INTEGER, INTENT (in) :: iq
    !! Current q-points
    INTEGER, INTENT (in) :: totq
    !! Total number of q-point from selecq
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
      IF ( iqq == 1 ) THEN
        OPEN(iufilfreq, file = filfreq, form = 'unformatted')
        WRITE(iufilfreq) totq, nmodes
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
    IF (iqq == 1) THEN
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
    ENDIF ! iqq
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
    IF ( iqq == 1 ) THEN 
       OPEN(iufileph, file = filephmat, form = 'unformatted')
    ELSE
       OPEN(iufileph, file = filephmat, position='append', form = 'unformatted')
    ENDIF
    !
    !IF ( iq .eq. 1 ) WRITE(iufileph,'(2i7)') my_pool_id+1, fermicount
    IF ( iqq == 1 ) WRITE(iufileph) my_pool_id+1, fermicount
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
    IF ( iqq == totq ) THEN 
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
       imelt = totq * nkfs
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
    SUBROUTINE count_kpoints( iqq )
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
    INTEGER, INTENT (in) :: iqq
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
    IF (iqq == 1) THEN
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
    ENDIF ! iqq
    !
    RETURN
    !
    END SUBROUTINE count_kpoints
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
    USE superconductivity, ONLY : mem_size_eliashberg, mem_integer_size_eliashberg
    USE constants_epw, ONLY : zero
    USE mp_global, ONLY : inter_pool_comm, npool
    USE mp,        ONLY : mp_bcast, mp_barrier, mp_sum
    USE mp_world,  ONLY : mpime
    USE division,  ONLY : fkbounds
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
    memlt_pool(:) = zero
#else
    IF ( .not. ALLOCATED(memlt_pool) ) ALLOCATE(memlt_pool(1))
    memlt_pool(1) = zero
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
    xkf_(:,:) = zero
    ekf_(:,:) = zero
    wkf_(:) = zero
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
    xkfs(:,:) = zero
    wkfs(:) = zero
    ekfs(:,:) = zero
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
    USE superconductivity, ONLY : mem_size_eliashberg, mem_integer_size_eliashberg
    USE constants_epw, ONLY : eps5, zero
    USE symm_base, ONLY : nrot
    USE io_global, ONLY : ionode_id
    USE mp_global, ONLY : inter_pool_comm
    USE mp,        ONLY : mp_bcast, mp_barrier, mp_sum
    USE mp_world,  ONLY : mpime
    USE division,  ONLY : fkbounds
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
    !! Lower bounds index after k paral
    INTEGER :: upper_bnd
    !! Upper bounds index after k paral
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
    xkff(:,:) = zero
    ixkff(:) = 0
    !
    ! to map k+q onto k we need to define the index of k on the full mesh (ixkff) 
    ! using index of the k-point within the Fermi shell (ixkf)
    !
    IF ( mpime .eq. ionode_id ) THEN
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
        CALL set_sym_bl( )
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
    USE constants_epw, ONLY : kelvin2eV, zero
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT (in) :: itemp
    !! Counter on temperature
    CHARACTER (len=256), INTENT (in) :: cname
    !
    ! Local variables
    INTEGER :: ik
    !! Counter on k-point
    INTEGER :: ibnd
    !! Counter on band
    INTEGER :: ibin
    !! Counter on bins
    INTEGER :: nbin
    !! Number of bins
    !
    REAL(DP) :: temp
    !! Temperature in K
    REAL(DP) :: dbin
    !! Step size in nbin
    REAL(DP) :: delta_max
    !! Max value of superconducting gap
    REAL(DP) :: sigma
    !! Variable for smearing
    REAL(DP) :: weight
    !! Variable for weight
    REAL(DP), ALLOCATABLE :: delta_k_bin(:)
    !! Histogram superconducting gap
    REAL(DP), EXTERNAL :: w0gauss
    CHARACTER (len=256) :: name1
    !
    temp = estemp(itemp) / kelvin2eV
    !
    delta_max = 1.25d0 * maxval(Agap(:,:,itemp))
    nbin = int(delta_max/(0.005d0/1000.d0))
    dbin = delta_max / dble(nbin)
    IF ( .not. ALLOCATED(delta_k_bin) ) ALLOCATE( delta_k_bin(nbin) )
    delta_k_bin(:) = zero
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
    USE constants_epw, ONLY : kelvin2eV, zero
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT (in) :: itemp
    !! Counter on temperature
    !
    ! Local variables
    INTEGER :: ik
    !! Counter on k-points
    INTEGER :: ibnd
    !! Counter on bands at k
    INTEGER :: i
    !! Counter on grid points nkf1
    INTEGER :: j
    !! Counter on grid points nkf2
    INTEGER :: k
    !! Counter on grid points nkf3
    !
    REAL(DP) :: temp
    !! Temperature in K
    REAL(DP) :: x1
    !! Cartesian coordinates of grid points nkf1
    REAL(DP) :: x2
    !! Cartesian coordinates of grid points nkf2
    REAL(DP) :: x3
    !! Cartesian coordinates of grid points nkf3
    REAL(DP), ALLOCATABLE :: Agap_tmp(:,:)
    !! Temporary array for superconducting gap at ik, ibnd
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
    Agap_tmp(:,0) = zero
    !
    ! SP & RM: .cube file for VESTA plotting (only if iverbosity = 2)
    !
    IF ( iverbosity .eq. 2 ) THEN
      !
      DO ibnd = 1, nbndfs
        !
        IF ( ibnd < 10 ) THEN
          ! We make the assumption that there are no superconductor with Tc
          ! higher than 999 K.
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
    ! SP & RM : Write on file the superconducting gap close to the Fermi surface
    ! along with
    !     Cartesian coordinate, band index, energy distance from Fermi level and
    !     gap value.
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
    !----------------------------------------------------------------------
    ! 
  END MODULE io_eliashberg


