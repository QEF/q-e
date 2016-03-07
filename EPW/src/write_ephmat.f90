  !                                                                            
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino 
  ! Copyright (C) 2007-2009 Roxana Margine 
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  !-----------------------------------------------------------------------
  SUBROUTINE write_ephmat( iq )
  !-----------------------------------------------------------------------
  !
  !  This subroutine writes the elph matrix elements in a format required 
  !  by Eliashberg equations
  ! 
  !  Use matrix elements, electronic eigenvalues and phonon frequencies
  !  from ep-wannier interpolation
  !
  !
  !-----------------------------------------------------------------------
#include "f_defs.h"
  USE kinds,      ONLY : DP
  USE io_global,  ONLY : stdout
  USE io_epw,     ONLY : iunepmatf, iufilfreq, iufilegnv, iufileph, iuetf
  USE io_files,   ONLY : prefix, tmp_dir
  USE phcom,      ONLY : nmodes
  USE control_lr, ONLY : lgamma
  USE epwcom,     ONLY : nbndsub, lrepmatf, fsthick, ngaussw, degaussw, & 
                         etf_mem, nkf1, nkf2, nkf3, eps_acustic, &
                         efermi_read, fermi_energy
  USE pwcom,      ONLY : nelec, ef, isk
  USE elph2,      ONLY : etf, ibndmin, ibndmax, nkqf, epf17, wkf, nkf, nqtotf, wf, xqf, nkqtotf
  USE eliashbergcom, ONLY : equivk, nkfs, ekfs, wkfs, xkfs, ef0, dosef, ixkf, ixkqf, nbndfs
  USE constants_epw, ONLY : ryd2ev, two
  USE mp_global,  ONLY :  my_pool_id
#ifdef __PARA
  USE mp,         ONLY : mp_barrier, mp_sum
  USE mp_global,  ONLY : inter_pool_comm, my_pool_id, npool
  USE mp_world,   ONLY : mpime
#endif
  !
  IMPLICIT NONE
  INTEGER :: ik, ikk, ikq, ibnd, jbnd, imode, iq, nrec, fermicount, & 
             nks, lower_bnd, upper_bnd, nkftot, imelt
  REAL(DP) :: wq, g2
  COMPLEX(kind=DP) epf (ibndmax-ibndmin+1, ibndmax-ibndmin+1)
  REAL(DP), EXTERNAL :: efermig, dos_ef
  CHARACTER (len=256) :: filfreq, filegnv, filephmat
#ifdef __PARA
  CHARACTER (len=3) :: filelab
#endif
  !
  ! write phonon frequencies to file
#ifdef __PARA 
  IF ( my_pool_id .eq. 0 ) THEN
#endif
  filfreq = trim(tmp_dir) // trim(prefix) // '.freq' 
  IF ( iq .eq. 1 ) THEN
     !OPEN(iufilfreq, file = filfreq, form = 'formatted')
     !WRITE(iufilfreq,'(2i7)') nqtotf, nmodes
     !WRITE(iufilfreq,'(3f15.9)') xqf(1,iq), xqf(2,iq), xqf(3,iq)
     !WRITE(iufilfreq,'(20ES20.10)') (wf(imode,iq), imode=1,nmodes)
     OPEN(iufilfreq, file = filfreq, form = 'unformatted')
     WRITE(iufilfreq) nqtotf, nmodes
     WRITE(iufilfreq) xqf(1,iq), xqf(2,iq), xqf(3,iq)
     DO imode = 1, nmodes
        WRITE(iufilfreq) wf(imode,iq)
     ENDDO
     CLOSE(iufilfreq)
  ELSE
     !OPEN(iufilfreq, file = filfreq, access = 'append', form = 'formatted')
     !WRITE(iufilfreq,'(3f15.9)') xqf(1,iq), xqf(2,iq), xqf(3,iq)
     !WRITE(iufilfreq,'(20ES20.10)') (wf(imode,iq), imode=1,nmodes)
     OPEN(iufilfreq, file = filfreq, position='append', form = 'unformatted')
     WRITE(iufilfreq) xqf(1,iq), xqf(2,iq), xqf(3,iq)
     DO imode = 1, nmodes
        WRITE(iufilfreq) wf(imode,iq)
     ENDDO
     CLOSE(iufilfreq)
  ENDIF
#ifdef __PARA
  ENDIF
#endif
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
        IF (.not. etf_mem) then
           nrec = ikk
           CALL davcio( etf(ibndmin:ibndmax,ikk), ibndmax-ibndmin+1, iuetf, nrec, -1 )
           nrec = ikq
           CALL davcio( etf(ibndmin:ibndmax,ikq), ibndmax-ibndmin+1, iuetf, nrec, -1 )
        ENDIF
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
#ifdef __PARA
     !
     ! collect contributions from all pools (sum over k-points)
     CALL mp_sum( nks, inter_pool_comm )
     CALL mp_barrier(inter_pool_comm)
     !
#endif
     !
     ! write eigenvalues to file
#ifdef __PARA 
  IF ( my_pool_id .eq. 0 ) THEN
#endif 
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
#ifdef __PARA
  ENDIF
#endif
     !
  ENDIF ! iq
  !
  ! write the e-ph matrix elements in the Bloch representation on the fine mesh
  ! in .ephmat files (one for each pool)
  !
#ifdef __PARA
  CALL set_ndnmbr(0,my_pool_id+1,1,npool,filelab)
  filephmat = trim(tmp_dir) // trim(prefix) // '.ephmat' // filelab
#else
  filephmat = trim(tmp_dir) // trim(prefix) // '.ephmat'
#endif
  IF ( iq .eq. 1 ) THEN 
     !OPEN(iufileph, file = filephmat, form = 'formatted')
     OPEN(iufileph, file = filephmat, form = 'unformatted')
  ELSE
     !OPEN(iufileph, file = filephmat, access = 'append', form = 'formatted')
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
     IF (lgamma) THEN
        ikk = ik
        ikq = ik
     ELSE
        ikk = 2 * ik - 1
        ikq = ikk + 1
     ENDIF
     !
     ! go only over irreducible k-points
     !
     IF ( equivk(lower_bnd+ik-1) .eq. (lower_bnd+ik-1) ) THEN 
        !
        ! we read the hamiltonian eigenvalues (those at k+q depend on q!)
        ! when we see references to iq for file reading, it is always = 1 
        IF (.not. etf_mem) then
           nrec = ikk
           CALL davcio( etf(ibndmin:ibndmax,ikk), ibndmax-ibndmin+1, iuetf, nrec, -1 )
           nrec = ikq
           CALL davcio( etf(ibndmin:ibndmax,ikq), ibndmax-ibndmin+1, iuetf, nrec, -1 )
        ENDIF
        !
        ! here we must have ef, not ef0, to be consistent with ephwann_shuffle
        !
        !   IF ( ixkf(lower_bnd+ik-1) .gt. 0 .AND. ixkqf(ixkf(lower_bnd+ik-1),iq) .gt. 0 ) THEN
        ! FG: here it can happen that ixkf is 0 and this leads to ixqf(0,iq) after .and.
        !     modified to prevent crash
        IF ( ixkf(lower_bnd+ik-1) .gt. 0 ) THEN
         IF ( ixkqf(ixkf(lower_bnd+ik-1),iq) .gt. 0 ) THEN

           DO imode = 1, nmodes ! phonon modes
              wq = wf(imode, iq)
              !
              !  we read the e-p matrix from memory / disk
              IF (etf_mem) THEN
                 epf(:,:) = epf17( ik, :, :, imode )
              ELSE
                 nrec = (imode-1) * nkf + ik
                 CALL dasmio( epf, ibndmax-ibndmin+1, lrepmatf, iunepmatf, nrec, -1 )
              ENDIF
              DO ibnd = 1, ibndmax-ibndmin+1
                 IF ( abs( ekfs(ibnd,ixkf(lower_bnd+ik-1)) - ef0 ) .lt. fsthick ) THEN
                    DO jbnd = 1, ibndmax-ibndmin+1
                       IF ( abs( ekfs(jbnd,ixkqf(ixkf(lower_bnd+ik-1),iq)) - ef0 ) .lt. fsthick ) THEN
                          !
                          ! here we take into account the zero-point sqrt(hbar/2M\omega)
                          ! with hbar = 1 and M already contained in the eigenmodes
                          ! g2 is Ry^2, wkf must already account for the spin factor
                          !
                          g2 = abs( epf(jbnd, ibnd) )**two / ( two * wq )
                          !WRITE(iufileph,'(ES20.10)') g2
                          WRITE(iufileph) g2
                       ENDIF
                    ENDDO ! jbnd
                 ENDIF
              ENDDO ! ibnd
           ENDDO ! imode

         ENDIF
        ENDIF ! fsthick

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
#include "f_defs.h"
  USE kinds,     ONLY : DP
  USE io_global, ONLY : stdout
  USE io_epw,    ONLY : iuetf
  USE epwcom,    ONLY : nbndsub, fsthick, ngaussw, degaussw, & 
                        etf_mem, efermi_read, fermi_energy, mp_mesh_k
  USE pwcom,     ONLY : nelec, ef, isk
  USE elph2,     ONLY : etf, ibndmin, ibndmax, nkqf, wkf, nkf, nkqtotf
  USE constants_epw, ONLY : two
#ifdef __PARA
  USE mp,        ONLY : mp_barrier, mp_sum
  USE mp_global, ONLY : inter_pool_comm, npool
  USE mp_world,  ONLY : mpime
#endif
  !
  IMPLICIT NONE
  INTEGER :: ik, ikk, ikq, iq, nrec, fermicount, nks
  REAL(DP) :: ef0, dosef
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
        IF (.not. etf_mem) then
           nrec = (iq-1) * nkqf + ikk
           nrec = ikk
           CALL davcio( etf(ibndmin:ibndmax,ikk), ibndmax-ibndmin+1, iuetf, nrec, -1 )
           nrec = (iq-1) * nkqf + ikq
           nrec = ikq
           CALL davcio( etf(ibndmin:ibndmax,ikq), ibndmax-ibndmin+1, iuetf, nrec, -1 )
        ENDIF
        !
        IF ( minval( abs( etf(:,ikk) - ef  ) ) .lt. fsthick ) &
           fermicount = fermicount + 1 
        !
     ENDDO
     !
     ! nks =  nr of k-points within the Fermi shell (fine mesh)
     nks = fermicount
     !
#ifdef __PARA
     !
     ! collect contributions from all pools (sum over k-points)
     CALL mp_sum( nks, inter_pool_comm )
     CALL mp_barrier(inter_pool_comm)
     !
#endif
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
