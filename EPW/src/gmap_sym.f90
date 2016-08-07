  !                                                                            
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino 
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  !-----------------------------------------------------------------------
  SUBROUTINE gmap_sym ( nsym, s, ftau, gmapsym, eigv, invs )
  !-----------------------------------------------------------------------
  !!
  !!   For every G vector, find S(G) for all the symmetry operations
  !!   of the crystal. Construct the matrix
  !!   eigv(ig,isym) = $e^{i G v(S)}$ where v(S) is the (possible) 
  !!   fractional translation associated with the symmetry operation
  !!
  !!   No parallelization on G-vecs at the moment  
  !!   (actually this is done on the global array, but in elphel2.f90
  !!   every processor has just a chunk of the array, I may need some
  !!   communication)
  !!
  !!   No ultrasoft now
  !!
  !----------------------------------------------------------------------
  USE kinds,         ONLY : DP
  USE pwcom,         ONLY : ngm
  USE constants_epw, ONLY : twopi, ci, cone
  USE fft_base,      ONLY : dfftp
  USE gvect,         ONLY : mill
  ! 
  implicit none
  !
  INTEGER, INTENT (in) :: nsym
  !! the number of symmetries of the crystal
  INTEGER, INTENT (in) :: s(3,3,48)
  !! the symmetry matrices
  INTEGER, INTENT (in) :: ftau(3,48)
  !! the fractional traslations
  INTEGER, INTENT (in) :: invs(48)
  !! inverse symmetry matrix 
  INTEGER, INTENT (out) :: gmapsym (ngm, 48)
  !! the map S(G) = gmapsym (G,S) 1...nsym
  COMPLEX(kind=DP), INTENT (out) :: eigv (ngm, 48)
  !! e^{ iGv} for 1...nsym
  !
  ! local variables
  !
  LOGICAL :: tfound
  !!
  INTEGER :: ig
  !! Counter on the G-vector
  INTEGER :: jg
  !! Counter on the G-vector
  INTEGER :: i
  !! Index of the rotated G-vector
  INTEGER :: j
  !! Index of the rotated G-vector
  INTEGER :: k
  !! Index of the rotated G-vector
  INTEGER :: notfound
  !! Index to check wether the mapping is complete
  INTEGER :: isym
  !! Counter on the symmetry 
  INTEGER :: ism1
  !! Index for the inverse symmetry
  REAL(kind=DP) :: rdotk
  !! $$\mathbf{r}\cdot\mathbf{k}
  !
  !  loop on the symmetries of the crystal
  !
  DO isym = 1, nsym
    !
    ism1 = invs(isym)
    !
    ! loop on the G vectors 
    !
    notfound = 0
    !
    DO ig = 1, ngm
      !
      !  the rotated G-vector
      !
      i = s (1, 1, isym) * mill(1,ig) + s (1, 2, isym) * mill(2,ig) + s (1, 3, isym) * mill(3,ig)
      j = s (2, 1, isym) * mill(1,ig) + s (2, 2, isym) * mill(2,ig) + s (2, 3, isym) * mill(3,ig)
      k = s (3, 1, isym) * mill(1,ig) + s (3, 2, isym) * mill(2,ig) + s (3, 3, isym) * mill(3,ig)
      !
      jg = 0
      tfound = .false.
      DO while ((.not.tfound).and.(jg.lt.ngm))
        jg = jg + 1
        tfound = (i.eq.mill(1,jg)) .and. (j.eq.mill(2,jg)) .and.(k.eq.mill(3,jg))
      ENDDO
      !
      IF (tfound) THEN
        gmapsym ( ig, isym ) = jg
      ELSE
        gmapsym ( ig, isym ) = 0
        notfound = notfound + 1
      ENDIF
      !
      ! now the phase factors e^{iGv}
      !
      IF ( ftau (1, isym).ne.0 .or. ftau (2, isym).ne.0 .or. ftau (3, isym).ne.0 ) THEN
        !
        ! fractional traslation in crystal coord is ftau/nr*
        ! for cart/crys transform of the G-vecctors have a look at the bottom
        !
        rdotk = dble( mill(1,ig) * ftau (1, isym) ) / dble (dfftp%nr1) &
              + dble( mill(2,ig) * ftau (2, isym) ) / dble (dfftp%nr2) &
              + dble( mill(3,ig) * ftau (3, isym) ) / dble (dfftp%nr3)
        !
        ! the actual translation is -v (have a look at ruota_ijk.f90)
        ! 
        eigv (ig, isym) = exp( - ci*twopi*rdotk ) 
        !
      ELSE
        eigv (ig, isym) = cone
      ENDIF
      !
    ENDDO
    !
    IF (notfound.gt.0) &
      CALL errore ('gmap_sym','incomplete mapping of G vectors: notfound = ',notfound)
    !
  ENDDO
  !
  END SUBROUTINE gmap_sym
