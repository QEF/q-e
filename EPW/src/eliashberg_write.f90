  !                                                                            
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino   
  ! Copyright (C) 2007-2009 Roxana Margine
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution_FS, or http://www.gnu.org/copyleft.gpl.txt .             
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
  USE constants_epw, ONLY : pi, kelvin2eV 
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
        WRITE(name1,'(a,a1,a4,a8,f4.2)') TRIM(prefix), '.', cname, '_aniso_0', temp
     ELSEIF ( temp .ge. 10.d0 ) THEN
        WRITE(name1,'(a,a1,a4,a7,f5.2)') TRIM(prefix), '.', cname, '_aniso_', temp
     ENDIF
     OPEN(iufilgap, file=name1, form='formatted')
     WRITE(iufilgap,'(5a20)') '#        w [eV]', 'Enk-Ef [eV]', 'Znorm(w) [eV]', 'Delta(w) [eV]', 'NZnorm(w) [eV]'
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
        WRITE(name1,'(a,a1,a4,a6,f4.2)') TRIM(prefix), '.', cname, '_iso_0', temp
     ELSEIF ( temp .ge. 10.d0 ) THEN
        WRITE(name1,'(a,a1,a4,a5,f5.2)') TRIM(prefix), '.', cname, '_iso_', temp
     ENDIF
     OPEN(iufilgap, file=name1, form='formatted')
     WRITE(iufilgap,'(4a24)') 'w', 'Znorm(w)', 'Delta(w)', 'NZnorm(w)'
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
           WRITE(name1,'(a,a1,a4,a8,f4.2)') TRIM(prefix), '.', cname, '_aniso_0', temp
        ELSEIF ( temp .ge. 10.d0 ) THEN
           WRITE(name1,'(a,a1,a4,a7,f5.2)') TRIM(prefix), '.', cname, '_aniso_', temp
        ENDIF
        OPEN(iufilgap, file=name1, form='formatted')
        WRITE(iufilgap,'(6a24)') 'w', 'Enk-Ef', 'Re[Znorm(w)]', 'Im[Znorm(w)]', 'Re[Delta(w)]', 'Im[Delta(w)]'
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
        WRITE(name1,'(a,a1,a4,a6,f4.2)') TRIM(prefix), '.', cname, '_iso_0', temp
     ELSEIF ( temp .ge. 10.d0 ) THEN
        WRITE(name1,'(a,a1,a4,a5,f5.2)') TRIM(prefix), '.', cname, '_iso_', temp
     ENDIF
     OPEN(iufilgap, file=name1, form='formatted')
     WRITE(iufilgap,'(5a18)') 'w', 'Re[Znorm(w)]', 'Im[Znorm(w)]', 'Re[Delta(w)]', 'Im[Delta(w)]'
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
     WRITE(name1,'(a,a1,a4,a13,f4.2)') TRIM(prefix), '.', cname, '_aniso_gap0_0', temp
  ELSEIF ( temp .ge. 10.d0 ) THEN
     WRITE(name1,'(a,a1,a4,a12,f5.2)') TRIM(prefix), '.', cname, '_aniso_gap0_', temp
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
  USE constants_epw, ONLY : pi, kelvin2eV
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
        IF ( temp .lt. 10.d0 ) THEN
           WRITE(name1,'(a,a1,a4,a13,f4.2,a1,i1,a5)')TRIM(prefix), '.', cname, '_aniso_gap0_0', temp, '_', ibnd, '.cube'
        ELSEIF ( temp .ge. 10.d0 ) THEN
           WRITE(name1,'(a,a1,a4,a12,f5.2,a1,i1,a5)')TRIM(prefix), '.', cname, '_aniso_gap0_', temp, '_', ibnd, '.cube'
        ENDIF
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
     WRITE(name1,'(a,a1,a4,a15,f4.2)') TRIM(prefix), '.', cname, '_aniso_gap_FS_0', temp
  ELSEIF ( temp .ge. 10.d0 ) THEN
     WRITE(name1,'(a,a1,a4,a14,f5.2)') TRIM(prefix), '.', cname, '_aniso_gap_FS_', temp
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
