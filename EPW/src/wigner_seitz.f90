  !                                                                            
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino 
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  ! Adapted from the original f77 Wannier code of Marzari, Vanderbilt, 
  !                    and Souza
  !
  !-----------------------------------------------------------------
  subroutine wigner_seitz (nk1, nk2, nk3, irvec, nrr, ndegen, wslen)
  !-----------------------------------------------------------------
  !!
  !! Calculates a grid of points that fall inside of (and eventually 
  !! on the surface of) the Wigner-Seitz supercell centered on the 
  !! origin of the Bravais lattice with primitive translations 
  !! nk1*a_1+nk2*a_2+nk3*a_3
  !!  
  !!
  !! w.r.t. the original version the wigner-seitz vectors are sorted
  !! by increasing lenght in output. In this way the electron and
  !! phonon wigner-seitz vectors should always be the same (even though
  !! the number of them may differ)
  !!
  !! BUG FIX: in the case of the tetragonal cell with c>a (LSCO)
  !! the WS points are correctly determined, but there are a few points
  !! with the wrong degeneracies. To avoid this I search for points
  !! in the -2:2 replicas (5^2 replicas). I had the same problem in 
  !! createkmap for the g0vec shift, and also there I have fixed it
  !! by extending the replicas to -2:2 instead of -1:1. FG May 07
  !!
  !-----------------------------------------------------------------
  USE kinds,     ONLY : DP
  USE cell_base, ONLY : at
  !
  implicit none
  !
  INTEGER, INTENT (in) :: nk1
  !! size of the uniform k mesh
  INTEGER, INTENT (in) :: nk2
  !! size of the uniform k mesh
  INTEGER, INTENT (in) :: nk3
  !! size of the uniform k mesh
  INTEGER, INTENT (out) :: irvec(3,20*nk1*nk2*nk3) 
  !! integer components of the ir-th Wigner-Seitz grid point in the basis of the lattice vectors
  INTEGER, INTENT (out) :: ndegen(20*nk1*nk2*nk3)
  !! Number of degeneracies
  INTEGER, INTENT (out) :: nrr
  !! number of Wigner-Seitz grid points 

  REAL(kind=DP), INTENT (out) :: wslen(20*nk1*nk2*nk3)
  !! real-space length, in units of alat
  !
  ! work variables
  integer :: irvec_ (3,20*nk1*nk2*nk3), ndegen_ (20*nk1*nk2*nk3)
  real(kind=DP), parameter :: eps = 1.d-8
  integer :: n1, n2, n3, i1, i2, i3, i, ipol, jpol, ndiff(3)!, ind(2*nk1*nk2*nk3)
  integer, allocatable :: ind(:)
  real(kind=DP) :: tot, mindist, adot(3,3), dist(125)
  logical :: found
  !
  !  the metric tensor 
  !
  INTEGER :: nind
  !
  nind = 20*nk1*nk2*nk3
  IF (nind .lt. 125) then
     nind = 125
  ENDIF
  IF (allocated(ind)) deallocate (ind)
  allocate (ind(nind))
  !
  DO ipol = 1, 3
   DO jpol = 1, 3
     adot (ipol, jpol) = dot_product ( at(:,ipol), at(:,jpol) )
   ENDDO
  ENDDO
  !
  ! Loop over grid points r on a unit cell that is 8 times larger than a 
  ! primitive supercell. In the end nrr contains the total number of grids 
  ! points that have been found in the Wigner-Seitz cell
  !
  nrr = 0
  DO n1 = 0, 4*nk1 
  DO n2 = 0, 4*nk2 
  DO n3 = 0, 4*nk3 
    !
    ! Loop over the 5^3 = 125 points R. R=0 corresponds to i1=i2=i3=2, or icnt=63
    !
    i = 0
    DO i1 = 0, 4
    DO i2 = 0, 4
    DO i3 = 0, 4
      i = i + 1
      !
      ! Calculate distance squared |r-R|^2 
      !
      ndiff(1) = n1-i1*nk1
      ndiff(2) = n2-i2*nk2
      ndiff(3) = n3-i3*nk3
      dist(i) = 0.d0
      DO ipol = 1, 3
       DO jpol = 1, 3
          dist(i) = dist(i) + dble(ndiff(ipol))*adot(ipol,jpol)*dble(ndiff(jpol))
       ENDDO
      ENDDO
      !  
    ENDDO 
    ENDDO 
    ENDDO
    !
    ! Sort the 125 vectors R by increasing value of |r-R|^2
    !
    ! NOTA BENE: hpsort really sorts the dist vector
    ! while the original subroutine by MVS did not. Therefore,
    ! dist(ind(i)) of the original version is here replacerd by
    ! dist(i), while ind(i) is kept.
    !
    ind(1) = 0 ! required for hpsort_eps (see the subroutine)
    CALL hpsort_eps_epw( 125, dist, ind, eps)

    !
    ! Find all the vectors R with the (same) smallest |r-R|^2;
    ! if R=0 is one of them, then the current point r belongs to 
    ! Wignez-Seitz cell => set found to true
    !
    found = .false.
    i = 1
    mindist = dist(1) 
    DO while ( abs(dist(i)-mindist).le.eps .and. i.lt.125 ) 
      IF (ind(i).eq.63) found = .true. 
      i = i + 1
    ENDDO
!@
    IF (i .eq. 126) i = 125
!@
    IF (found) then
      nrr = nrr + 1
      ndegen (nrr) = i - 1
      irvec (1, nrr) = n1 - 2*nk1
      irvec (2, nrr) = n2 - 2*nk2
      irvec (3, nrr) = n3 - 2*nk3
    ENDIF
    !
  ENDDO   
  ENDDO  
  ENDDO 
  !
  ! Check the "sum rule"
  !
  tot = 0.d0
  DO i = 1, nrr
   tot = tot + 1.d0 / dble (ndegen(i)) 
  ENDDO
  !
  IF(abs(tot-dble(nk1*nk2*nk3)).gt.eps) call errore &
       ('wigner_seitz','weights do not add up to nk1*nk2*nk3',1)
  !
  !@ JN it happens in 2d and 1d systems with small course grids.  I've changed to 20**3
  !@ could calculate the max number of elements at the beginning
  ! Hopefully this will never happen, i.e., I think 2*nk1*nk2*nk3 is
  ! an upper bound to the number of lattice points in (or on
  ! the surface of) the Wigner-Seitz supercell
  !
  IF(nrr.gt.20*nk1*nk2*nk3) call errore &
    ('wigner_seitz','too many wigseit points, try to increase the bound 20*nk1*nk2*nk3',1)
  !
  ! Now sort the wigner-seitz vectors by increasing magnitude
  !
  DO i = 1, nrr
    wslen(i) = 0.d0
    DO ipol = 1, 3
     DO jpol = 1, 3
        wslen(i) = wslen(i) + dble(irvec(ipol,i))*adot(ipol,jpol)*dble(irvec(jpol,i))
     ENDDO
    ENDDO
    wslen(i) = sqrt(wslen(i))
  ENDDO 
  !
  ind(1) = 0 ! required for hpsort_eps (see the subroutine)
  CALL hpsort_eps( nrr, wslen, ind, eps)
  !
  !  now wslen is already sorted, but we still have to sort
  !  irvec and ndegen
  !
  DO i = 1, nrr
    ndegen_ (i)  = ndegen(ind(i))
    irvec_ (:,i) = irvec(:,ind(i))
  ENDDO
  ndegen = ndegen_
  irvec  = irvec_
  !
  end subroutine wigner_seitz

