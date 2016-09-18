  !                                                                            
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino 
  ! Copyright (C) 2007-2009 Roxana Margine
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
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
  USE io_global, ONLY : stdout
  USE io_epw,    ONLY : iufilfreq
  USE io_files,  ONLY : prefix, tmp_dir
  USE phcom,     ONLY : nmodes
  USE elph2,   ONLY : nqtotf, wf, wqf, xqf
  USE eliashbergcom, ONLY : wsphmax
  USE constants_epw, ONLY : ryd2ev
  USE io_global, ONLY : ionode_id 
  USE mp_global, ONLY : inter_pool_comm
  USE mp_world,  ONLY : mpime
  USE mp,        ONLY : mp_bcast, mp_barrier, mp_sum
  !
  IMPLICIT NONE
  !    
  INTEGER :: iq, imode, ios
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
  USE io_global,     ONLY : stdout
  USE io_files,      ONLY : prefix, tmp_dir
  USE pwcom,         ONLY : ef
  USE epwcom,        ONLY : nkf1, nkf2, nkf3, degaussw, fsthick, mp_mesh_k
  USE eliashbergcom, ONLY : nkfs, nbndfs, dosef, ef0, ekfs, wkfs, xkfs, w0g
  USE constants_epw, ONLY : ryd2ev
  USE io_epw,        ONLY : iufilegnv
  USE io_global, ONLY : ionode_id
  USE mp_global, ONLY : inter_pool_comm
  USE mp_world,  ONLY : mpime
  USE mp,        ONLY : mp_bcast, mp_barrier, mp_sum
  ! 
  IMPLICIT NONE
  !
  REAL(DP), ALLOCATABLE :: ekf_(:,:)
  INTEGER :: ik, ibnd, ios, nkftot, n, nbnd_
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
  USE io_global, ONLY : stdout
  USE io_epw,    ONLY : iufilikmap
  USE io_files,  ONLY : prefix, tmp_dir
  USE symm_base, ONLY : t_rev, time_reversal, s, set_sym_bl
  USE phcom,     ONLY : nmodes
  USE epwcom,    ONLY : nkf1, nkf2, nkf3, mp_mesh_k
  USE elph2,     ONLY : nqtotf, xqf
  USE eliashbergcom, ONLY : ixkff, xkff, ixkf, xkfs, nkfs, ixkqf, ixqfs, nbndfs, nqfs, memlt_pool
  USE mp_global,     ONLY : npool
  USE symm_base, ONLY : nrot
  USE io_global, ONLY : ionode_id
  USE mp_global, ONLY : inter_pool_comm
  USE mp_world,  ONLY : mpime
  USE mp,        ONLY : mp_bcast, mp_barrier, mp_sum
  ! 
  IMPLICIT NONE
  !
  REAL(DP), PARAMETER :: eps=1.0d-5
  INTEGER :: i, j, k, ik, iq, n, nkq, nks, nk, nkf_irr, nkftot, ns,& 
             lower_bnd, upper_bnd, ios, imelt
  INTEGER, ALLOCATABLE :: index_(:,:), equiv_(:)
  REAL(DP) :: xk(3), xq(3), xkr(3), xx, yy, zz
  LOGICAL :: in_the_list
  CHARACTER (len=256) :: filikmap
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
    ! nkf_irr - total nr of irreducible k-points
    !READ(iufilikmap,'(i9)') nkf_irr
    READ(iufilikmap) nkf_irr
    !
    IF ( .not. ALLOCATED(ixkf) ) ALLOCATE(ixkf(nkf_irr))
    ixkf(:) = 0
    !
    DO ik = 1, nkf_irr
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
       equiv_(nk)=nk
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
                in_the_list = abs( xx-nint(xx) ) .le. eps .AND. &
                              abs( yy-nint(yy) ) .le. eps .AND. &
                              abs( zz-nint(zz) ) .le. eps
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
                   in_the_list = abs( xx-nint(xx) ) .le. eps .AND. &
                                 abs( yy-nint(yy) ) .le. eps .AND. &
                                 abs( zz-nint(zz) ) .le. eps
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
    IF ( nks .ne. nkf_irr) CALL errore('read_kmap_mp', 'something wrong with the mesh',1)
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
  ! find the index of k+sign*q on the irreducible fine k-mesh
  ! nkfs - nr of irreducible k-points within the Fermi shell (fine mesh)
  ! nqtotf - total nr of q-points on the fine mesh
  !
  DO ik = lower_bnd, upper_bnd
     DO iq = 1, nqtotf
        xk(:) = xkfs(:,ik)
        xq(:) = xqf(:,iq)
        !
        !  nkq - index of k+sign*q on the irreducible fine k-mesh.
        !
        CALL kpmq_map( xk, xq, +1, nkq )
        !
        !  ixkqf(ik,iq) - index of k+sign*q on the irreducible fine k-mesh
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
  USE mp_global,     ONLY : npool
  USE mp,            ONLY : mp_barrier, mp_bcast, mp_sum
  USE mp_global,     ONLY : inter_pool_comm
  !  
  IMPLICIT NONE
  !
  INTEGER :: ik, iq, ibnd, jbnd, imode, nnk, nnq(nkfs), ipool, tmp_pool_id, ios, & 
             lower_bnd, upper_bnd, nkpool(npool), nmin, nmax, nks, imelt
  REAL(DP) :: gmat
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
