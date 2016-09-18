  !                                                                            
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino 
  ! Copyright (C) 2007-2009 Roxana Margine
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  !-----------------------------------------------------------------------
  SUBROUTINE kmesh_fine
  !-----------------------------------------------------------------------
  !!
  !!   This routine defines the nr. of k-points on the fine k-mesh 
  !!   within the Fermi shell
  !!
  USE kinds,     ONLY : DP
  USE io_global, ONLY : stdout
  USE io_files,  ONLY : prefix, tmp_dir
  USE epwcom,    ONLY : nkf1, nkf2, nkf3, fsthick, mp_mesh_k
  USE pwcom,     ONLY : ef
  USE io_epw,    ONLY : iufilikmap
  USE elph2,     ONLY : xkf, wkf, etf, nkf, nkqtotf, ibndmin, ibndmax
  USE eliashbergcom, ONLY : nkfs, ixkf, equivk, xkfs, wkfs, ekfs, nbndfs, memlt_pool
  USE io_global, ONLY : ionode_id
  USE mp_global, ONLY : inter_pool_comm, npool
  USE mp,        ONLY : mp_bcast, mp_barrier, mp_sum
  USE mp_world,  ONLY : mpime
  !
  IMPLICIT NONE
  !
  INTEGER :: nk, nks, ikk, lower_bnd, upper_bnd, nkf_mesh, imelt
  REAL(DP) :: xx, yy, zz
  REAL(DP), ALLOCATABLE :: wkf_(:), ekf_(:,:), xkf_(:,:)
  CHARACTER (len=256) :: filikmap
  !
  ! nkf_mesh - nr of k-points 
  ! for mp_mesh_k = true, nkf_mesh - nr of irreducible k-points
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
       equivk(nk)=nk
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
  USE symm_base, ONLY : nrot
  USE io_global, ONLY : ionode_id
  USE mp_global, ONLY : inter_pool_comm
  USE mp,        ONLY : mp_bcast, mp_barrier, mp_sum
  USE mp_world,  ONLY : mpime
  ! 
  IMPLICIT NONE
  !
  REAL(DP), PARAMETER :: eps=1.0d-5
  INTEGER :: i, j, k, ik, iq, n, nkq, nks, nk, nkftot, ns, lower_bnd, upper_bnd, imelt
  INTEGER, ALLOCATABLE :: index_(:,:), equiv_(:)
  REAL(DP) :: xk(3), xq(3), xkr(3), xx, yy, zz
  LOGICAL :: in_the_list
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
  ! nkfs - nr of k-points within the Fermi shell  
  ! for mp_mesh_k = true, nkfs - nr of irreducible k-points within the Fermi shell
  ! nqtotf - total nr of q-points on the fine mesh
  !
  DO ik = lower_bnd, upper_bnd
     DO iq = 1, nqtotf
        xk(:) = xkfs(:,ik)
        xq(:) = xqf(:,iq)
        !
        ! find nkq - index of k+sign*q on the irreducible fine k-mesh.
        !
        CALL kpmq_map( xk, xq, +1, nkq )
        !
        ! find ixkqf(ik,iq) - index of k+sign*q on the irreducible fine k-mesh
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
  ENDDO ! loop over irreducible k-points within the Fermi shell in each pool (fine mesh) 
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
  REAL(DP) :: xx, yy, zz, eps, xxk(3)
  LOGICAL :: in_the_list
  !
  !
  ! loosy tolerance, no problem since we use integer comparisons
  eps = 1.d-5
  !
  xxk(:) = xk(:) + dble(sign1) * xq(:)
  xx = xxk(1) * nkf1
  yy = xxk(2) * nkf2
  zz = xxk(3) * nkf3
  in_the_list = abs(xx-nint(xx)) .le. eps .AND. &
                abs(yy-nint(yy)) .le. eps .AND. &
                abs(zz-nint(zz)) .le. eps
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
