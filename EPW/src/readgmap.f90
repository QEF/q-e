  !                                                                            
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino 
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  !--------------------------------------------------------------
  SUBROUTINE readgmap( nkstot, ngxx, ng0vec, g0vec_all_r, lower_bnd ) 
  !--------------------------------------------------------------
  !!
  !!  read map of G vectors G -> G-G_0 for a given q point
  !!  (this is used for the folding of k+q into the first BZ) 
  !!    
  !!
  !--------------------------------------------------------------
  USE kinds,    ONLY : DP
  USE mp_global,ONLY : inter_pool_comm, world_comm
  USE mp,       ONLY : mp_bcast, mp_max
  use io_global,ONLY : meta_ionode, meta_ionode_id
  use io_epw,   ONLY : iukgmap, iukmap
  use pwcom,    ONLY : nks
  use elph2,    ONLY : shift, gmap, igk_k_all, ngk_all
  USE io_files, ONLY : prefix
  !
  IMPLICIT NONE
  !
  ! variables for folding of k+q grid
  !
  INTEGER, INTENT(in) :: nkstot
  !! Total number of k-points
  INTEGER, INTENT(out) :: ngxx
  !! Maximum number of G-vectors over all pools
  INTEGER, INTENT(out) :: ng0vec
  !! Number of G_0 vectors
  INTEGER, INTENT(in) :: lower_bnd
  !! Lower bound for the k-parallellization
  ! 
  REAL(kind=DP), INTENT(out) :: g0vec_all_r(3,125)
  !! G_0 vectors needed to fold the k+q grid into the k grid, cartesian coord.
  !
  !  work variables
  !
  INTEGER :: ik
  !! Counter on k-points 
  INTEGER :: ik1, itmp
  !! Temporary indeces when reading kmap and kgmap files
  INTEGER :: ig0
  !! Counter on G_0 vectors
  INTEGER :: ishift
  !! Counter on G_0 vectors
  INTEGER :: ig
  !! Counter on G vectors
  INTEGER :: ios
  !! Integer variable for I/O control
  !
  REAL(DP) :: tmp
  !
  ALLOCATE( shift(nkstot) )
  !
  !  OBSOLETE: now we read directly the igkq to get the proper ngxx
  !
  !  read only a piece of the map to save time 
  !  the proper allocation bound would be ngxx = max(max(igkq))
  !  where the max is taken over the ig and the ik
  !  Here I use a simpler estimate: take the sphere npwx + two
  !  extra shells. This may not work for strange shapes of the
  !  reciproc latt. In this case just set ngxx = ngm_g
  !
  !  ngxx = nint(4./3.*3.14*(2+(3.0/4.0/3.14*dble(npwx))**(1./3.))**3.)
  !

  !  Note that the k+q point below does not correspond to the actual (true) 
  !  k+q, but since we only need to take the max over k and k+q this
  !  does not matter
  !
  ngxx = 0
  DO ik = 1, nks
    !
    IF ( maxval(igk_k_all(1:ngk_all(ik+lower_bnd-1),ik+lower_bnd-1)) > ngxx ) &
      ngxx = maxval(igk_k_all(1:ngk_all(ik+lower_bnd-1),ik+lower_bnd-1))
    !
  ENDDO
  !
#if defined(__MPI)
  tmp = dble(ngxx)
  CALL mp_max( tmp, inter_pool_comm )  
  ngxx = nint(tmp)
#endif
  !
  IF (meta_ionode) THEN
    !
    OPEN(iukgmap, file=trim(prefix)//'.kgmap', form='formatted', status='old', iostat=ios)
    IF (ios /=0) CALL errore('readgmap', 'error opening kgmap file', iukgmap)
    !
    DO ik = 1, nkstot
      READ(iukgmap,*) ik1, shift(ik1)
    ENDDO
    READ(iukgmap,*) ng0vec
    !
    !  the following seems crazy but I make it for compatibility
    !  with versions up to 2.1.5:
    !
    !  iukgmap has been created by ../PW/set_kplusq.f90 and has
    !  the correct gmap(), but the wrong shift() (actually the
    !  shift for a specific q-point)
    !
    !  since createkmap.f90 has regenerated the shifts for the
    !  present k-point I read them again in kmap.dat. The above 
    !  'fake' reading is because the gmap appears *after* the
    !  wrong kmap.
    !
    OPEN(iukmap, file=trim(prefix)//'.kmap', form='formatted', status='old', iostat=ios)
    IF (ios /= 0) CALL errore ('readgmap', 'error opening kmap file', iukmap)
    DO ik = 1, nkstot
      READ(iukmap,*) ik1, itmp, shift(ik1)
    ENDDO
    CLOSE(iukmap) 
    !
  ENDIF
  !
  ! first node broadcasts ng0vec to all nodes for allocation of gmap
  !
  CALL mp_bcast( ng0vec, meta_ionode_id, world_comm )
  !
  ALLOCATE( gmap(ngxx * ng0vec) )
  !
  IF (meta_ionode) THEN
     !
    DO ig0 = 1, ng0vec
      READ(iukgmap,*) g0vec_all_r(:,ig0)
    ENDDO
    DO ig = 1, ngxx
      ! 
      ! at variance with the nscf calculation, here gmap is read as a vector,
      ! 
      READ(iukgmap,*) ( gmap(ng0vec * ( ig - 1 ) + ishift), ishift = 1, ng0vec )
    ENDDO
    !
    CLOSE(iukgmap)
    !
  ENDIF
  !
  ! first node broadcasts everything to all nodes
  !
  CALL mp_bcast( g0vec_all_r, meta_ionode_id, world_comm )
  CALL mp_bcast( shift, meta_ionode_id, world_comm )
  CALL mp_bcast( gmap, meta_ionode_id, world_comm )
  !
  END SUBROUTINE readgmap
