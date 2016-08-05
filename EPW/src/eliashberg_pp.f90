  !                                                                            
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino 
  ! Copyright (C) 2007-2009 Roxana Margine
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
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
     WRITE(fildos,'(a,a7,f4.2)') TRIM(prefix),'.qdos_0', temp
  ELSEIF ( temp .ge. 10.d0 ) THEN
     WRITE(fildos,'(a,a6,f5.2)') TRIM(prefix),'.qdos_', temp
  ENDIF
  OPEN(iuqdos, file=fildos, form='formatted')
  !
  IF ( .not. ALLOCATED(dos_qp) ) ALLOCATE( dos_qp(nsw) )
  dos_qp(:) = 0.d0          
  !
  IF ( laniso ) THEN
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
     WRITE(filfe,'(a,a5,f4.2)') TRIM(prefix),'.fe_0', temp
  ELSEIF ( temp .ge. 10.d0 ) THEN
     WRITE(filfe,'(a,a4,f5.2)') TRIM(prefix),'.fe_', temp
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

