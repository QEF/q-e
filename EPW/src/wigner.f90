  !                                                                            
  ! Copyright (C) 2010-2019 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino 
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  !----------------------------------------------------------------------
  MODULE wigner
  !----------------------------------------------------------------------
  !! 
  !! This module contains all the subroutine linked creation of Wigner-Seitz cell
  !! 
  IMPLICIT NONE
  ! 
  CONTAINS
    !-----------------------------------------------------------------
    SUBROUTINE wigner_seitz_wrap (nk1, nk2, nk3, nq1, nq2, nq3, &
                                  irvec_k,  irvec_q,  irvec_g,  &
                                  ndegen_k, ndegen_q, ndegen_g, &
                                  wslen_k,  wslen_q,  wslen_g, w_centers, dims)
    !-----------------------------------------------------------------
    !!
    !! June 2018 - SP - CV
    !! 
    !! This routine wrap the call to three Wigner-Seitz routines:
    !!   wigner_seitzk : Contruct a grid of points that fall inside of (and eventually 
    !!                   on the surface of) the Wigner-Seitz supercell centered on the 
    !!                   origin of the Bravais lattice. Use for electronic properties. 
    !!   wigner_seitzq : Creates a set of WS vectors for each pair of atoms tau(nb)-tau(na)
    !!                   On exiting, ndegen_q contains the degeneracies of each pairs 
    !!                   of atoms while irvec_q contains the minimal communal sets of WS vectors. 
    !!                   Used for phonon properties
    !!   wigner_seitzg : Creates a set of WS vector for each atoms tau(na). 
    !!                   On exiting, ndegen_g contains the degeneracies of each atoms 
    !!                   while irvec_g contains the minimal communal sets of WS vectors. 
    !!                   Used for electron-phonon properties.  
    !! 
    !! Note 1: ndegen_k is always > 0 while ndegen_q and ndegen_g might contains 0 weigths. 
    !! Note 2: No sorting of vectors is needed anymore  
    !! Note 3: The dimension 20*nk1*nk2*nk3 should be safe enough.
    !! Note 4: The Wigner-Seitz construction in EPW was done by constructing a cell
    !!         centred unit cell. This is fine for electronic properties (this is what is done in wannier90).
    !!         However for phonon or electron-phonon properties, one can have issues when the cell
    !!         is tilted for example.
    !!         The proper way is to construct a set of WS vectors centred on pairs of atoms (phonons)
    !!         or atoms (el-ph).
    !!         In the matdyn code, a FT grid is constructed with weights centred on pairs of atoms 
    !!         and zeros everywhere else.
    !!         EPW now reproduced exactly the results of matdyn for the interpolated phonons at a
    !!         lower computation cost. Indeed we minimize the number of zeros by keeping the union
    !!         of values between all the cells.
    !!         In both cases this is very fast anyway but is important for el-ph properties.
    !!
    !-----------------------------------------------------------------
    USE kinds,     ONLY : DP
    USE ions_base, ONLY : nat, tau
    !
    implicit none
    !
    INTEGER, INTENT (in) :: nk1
    !! size of the uniform k mesh
    INTEGER, INTENT (in) :: nk2
    !! size of the uniform k mesh
    INTEGER, INTENT (in) :: nk3
    !! size of the uniform k mesh
    INTEGER, INTENT (in) :: nq1
    !! size of the uniform q mesh
    INTEGER, INTENT (in) :: nq2
    !! size of the uniform q mesh
    INTEGER, INTENT (in) :: nq3
    !! size of the uniform q mesh
    INTEGER, INTENT (in) :: dims
    !! Number of bands in the Wannier space
    INTEGER, ALLOCATABLE, INTENT (out) :: irvec_k(:,:)
    !! integer components of the ir-th Wigner-Seitz grid point in the basis
    !! of the lattice vectors for electrons
    INTEGER, ALLOCATABLE, INTENT (out) :: irvec_q(:,:)
    !! integer components of the ir-th Wigner-Seitz grid point for phonons
    INTEGER, ALLOCATABLE, INTENT (out) :: irvec_g(:,:)
    !! integer components of the ir-th Wigner-Seitz grid point for electron-phonon
    INTEGER, ALLOCATABLE, INTENT (out) :: ndegen_k (:,:,:)
    !! Wigner-Seitz number of degenerescence (weights) for the electrons grid
    INTEGER, ALLOCATABLE, INTENT (out) :: ndegen_q (:,:,:)
    !! Wigner-Seitz weights for the phonon grid that depend on 
    !! atomic positions $R + \tau(nb) - \tau(na)$
    INTEGER, ALLOCATABLE, INTENT (out) :: ndegen_g (:,:,:,:)
    !! Wigner-Seitz weights for the electron-phonon grid that depend on 
    !! atomic positions $R - \tau(na)$
    REAL(kind=DP), ALLOCATABLE, INTENT (out) :: wslen_k(:)
    !! real-space length for electrons, in units of alat
    REAL(kind=DP), ALLOCATABLE, INTENT (out) :: wslen_q(:)
    !! real-space length for phonons, in units of alat
    REAL(kind=DP), ALLOCATABLE, INTENT (out) :: wslen_g(:)
    !! real-space length for electron-phonons, in units of alat
    REAL(kind=DP), INTENT(in) :: w_centers(3,dims)
    !! Wannier centers used for the creation of electron and el-ph WS cells
    ! 
    ! Work Variables
    INTEGER :: ir
    !! Index for WS vectors
    INTEGER :: nrr_k
    !! maximum number of WS vectors for the electrons
    INTEGER :: nrr_q
    !! maximum number of WS vectors for the phonons
    INTEGER :: nrr_g
    !! maximum number of WS vectors for the electron-phonon
    INTEGER :: irvec_kk (3,20*nk1*nk2*nk3)
    !! local integer components of the ir-th Wigner-Seitz grid point 
    !! in the basis of the lattice vectors for electrons
    INTEGER :: irvec_qq (3,20*nq1*nq2*nq3)
    !! local integer components of the ir-th Wigner-Seitz grid point for phonons
    INTEGER :: irvec_gg (3,20*nq1*nq2*nq3)
    !! local integer components of the ir-th Wigner-Seitz grid point for electron-phonons
    !! We use nk1 instead of nq1 because the k-grid is always larger or equal to q-grid.  
    INTEGER :: ndegen_kk (20*nk1*nk2*nk3, dims, dims)
    !! local Wigner-Seitz number of degenerescence (weights) for the electrons grid
    INTEGER :: ndegen_qq (20*nq1*nq2*nq3, nat, nat)
    !! local Wigner-Seitz number of degenerescence (weights) for the phonons grid
    INTEGER :: ndegen_gg (20*nq1*nq2*nq3, nat, dims, dims)
    !! local Wigner-Seitz number of degenerescence (weights) for the electron-phonons grid
    REAL(kind=DP) :: wslen_kk (20*nk1*nk2*nk3)
    !! local real-space length for electrons, in units of alat
    REAL(kind=DP) :: wslen_qq (20*nq1*nq2*nq3)
    !! local real-space length for phonons, in units of alat
    REAL(kind=DP) :: wslen_gg (20*nq1*nq2*nq3)
    !! local real-space length for electron-phonon, in units of alat
    !
    !  Check the bounds
    IF ( nq1 > nk1 .OR. nq2 > nk2 .OR. nq3 > nk3 ) call errore &
       ('wigner_seitz_wrap',' the phonon grid should be smaller than electron grid',1)
    !
    ! Now generated the un-sorted points for the electrons, phonons and electron-phonon
    !
    ! If dims > 1, it includes the position of Wannier-Centers
    CALL wigner_seitzkq ( nk1, nk2, nk3, irvec_kk, ndegen_kk, wslen_kk, nrr_k, w_centers, dims)
    CALL wigner_seitzkq ( nq1, nq2, nq3, irvec_qq, ndegen_qq, wslen_qq, nrr_q, tau, nat)
    CALL wigner_seitzg  ( nq1, nq2, nq3, irvec_gg, ndegen_gg, wslen_gg, nrr_g, w_centers, tau, dims, nat)
    ! 
    ALLOCATE ( irvec_k(3,nrr_k) )
    ALLOCATE ( irvec_q(3,nrr_q) )
    ALLOCATE ( irvec_g(3,nrr_g) )
    ALLOCATE ( ndegen_k(nrr_k, dims, dims) )
    ALLOCATE ( ndegen_q(nrr_q, nat, nat) )
    ALLOCATE ( ndegen_g(nrr_g, nat, dims, dims) )
    ALLOCATE ( wslen_k(nrr_k) )
    ALLOCATE ( wslen_q(nrr_q) )
    ALLOCATE ( wslen_g(nrr_g) )
    ! 
    ! Create vectors with correct size. 
    DO ir = 1, nrr_k
      ndegen_k(ir,:,:)  = ndegen_kk(ir,:,:)
      irvec_k(:,ir)     = irvec_kk(:,ir)
      wslen_k(ir)       = wslen_kk(ir)
    ENDDO
    DO ir = 1, nrr_q
      ndegen_q(ir,:,:)  = ndegen_qq(ir,:,:)
      irvec_q(:,ir)     = irvec_qq(:,ir)
      wslen_q(ir)       = wslen_qq(ir)  
    ENDDO
    DO ir = 1, nrr_g
      ndegen_g(ir,:,:,:)  = ndegen_gg(ir,:,:,:)
      irvec_g(:,ir)       = irvec_gg(:,ir)
      wslen_g(ir)         = wslen_gg(ir)  
    ENDDO
    !
    END SUBROUTINE wigner_seitz_wrap
    !-----------------------------------------------------------------------------
    ! 
    ! 
    !-----------------------------------------------------------------------------
    SUBROUTINE wigner_seitzkq (nc1, nc2, nc3, irvec, ndegen, wslen, nrr, shift, dims)
    !-----------------------------------------------------------------------------    
    !!
    !! Calculates a grid of points that fall inside of (and eventually 
    !! on the surface of) the Wigner-Seitz supercell centered on the 
    !! origin of the Bravais lattice with primitive translations 
    !! nk1*a_1+nk2*a_2+nk3*a_3
    !!  
    !-----------------------------------------------------------------
    USE kinds,         ONLY : DP
    USE cell_base,     ONLY : at, bg
    USE constants_epw, ONLY : eps8
    !
    implicit none
    !
    INTEGER, INTENT (in) :: nc1
    !! size of the uniform k mesh
    INTEGER, INTENT (in) :: nc2
    !! size of the uniform k mesh
    INTEGER, INTENT (in) :: nc3
    !! size of the uniform k mesh
    INTEGER, INTENT (in) :: dims
    !! dim of shift(3,dims)
    INTEGER, INTENT (out) :: irvec(3,20*nc1*nc2*nc3)
    !! integer components of the ir-th Wigner-Seitz grid point in the basis of the lattice vectors
    INTEGER, INTENT (out) :: ndegen(20*nc1*nc2*nc3, dims, dims)
    !! Number of degeneracies
    INTEGER, INTENT (out) :: nrr
    !! number of Wigner-Seitz grid points 
    !
    REAL(kind=DP), INTENT (out) :: wslen(20*nc1*nc2*nc3)
    !! real-space length, in units of alat
    REAL(kind=DP), INTENT(in) :: shift(3,dims)
    !! Wannier centers (electron WS cell) or atomic positions (lattice WS cell)
    !
    ! Work variables
    INTEGER :: n1, n2, n3
    !! Index for the larger supercell
    INTEGER :: i1, i2, i3
    !! Index to compute |r-R| distance
    INTEGER :: i, ir, irtot, iw, iw2
    !! Iterative index
    INTEGER :: ipol, jpol
    !! Cartesian direction 
    INTEGER, ALLOCATABLE :: ind(:)
    !! Index of sorting
    INTEGER :: nind 
    !! The metric tensor
    INTEGER :: nrr_tmp(dims,dims)
    !! Temporary WS matrix
    INTEGER :: irvec_tmp(3,20*nc1*nc2*nc3,dims,dims)
    !! Temp
    INTEGER :: ndegen_tmp(20*nc1*nc2*nc3,dims,dims) 
    ! 
    REAL(kind=DP) :: adot(3,3) 
    !! Dot product between lattice vector
    REAL(kind=DP) :: dist(125)
    !! Contains the distance squared |r-R|^2
    REAL(kind=DP) :: mindist
    !! Minimum distance
    REAL(kind=DP) :: tot, tot2
    !! Sum of all the weigths
    REAL(kind=DP) :: ndiff(3)
    !! Distances. Must be now real because of fractional atoms 
    ! 
    LOGICAL :: found
    !! True if the vector has been found
    !
    nind = 20*nc1*nc2*nc3
    IF (nind < 125) then
      nind = 125
    ENDIF
    ALLOCATE (ind(nind))
    ! 
    DO ipol = 1, 3
     DO jpol = 1, 3
       adot (ipol, jpol) = dot_product ( at(:,ipol), at(:,jpol) )
     ENDDO
    ENDDO
    !
    CALL cryst_to_cart(dims,shift(:,:),bg,-1)
    !
    ndegen_tmp(:,:,:) = 0
    irvec_tmp(:,:,:,:) = 0
    DO iw = 1, dims
      DO iw2 = 1, dims
        nrr_tmp(iw,iw2) = 0
        DO n1 = -2*nc1, 2*nc1
          DO n2 = -2*nc2, 2*nc2
            DO n3 = -2*nc3, 2*nc3
              ! Loop over the 5^3 = 125 points R. R=0 corresponds to i1=i2=i3=2, or icnt=63
              i = 0
              dist(:) = 0.d0
              DO i1 = -2, 2
                DO i2 = -2, 2
                  DO i3 = -2, 2
                    i = i + 1  
                    ! Calculate distance squared |r-R|^2
                    ndiff(1) = n1 - i1*nc1 + shift(1,iw2) - shift(1,iw)
                    ndiff(2) = n2 - i2*nc2 + shift(2,iw2) - shift(2,iw)
                    ndiff(3) = n3 - i3*nc3 + shift(3,iw2) - shift(3,iw)
                    !ndiff(1) = n1 - i1*nc1 + (shift(1,iw2) + shift(1,iw)) / 2.d0
                    !ndiff(2) = n2 - i2*nc2 + (shift(2,iw2) + shift(2,iw)) / 2.d0
                    !ndiff(3) = n3 - i3*nc3 + (shift(3,iw2) + shift(3,iw)) / 2.d0
                    DO ipol = 1, 3
                      DO jpol = 1, 3
                        dist(i) = dist(i) + dble(ndiff(ipol))*adot(ipol,jpol)*dble(ndiff(jpol))
                      ENDDO
                    ENDDO
                    !  
                  ENDDO ! i3
                ENDDO ! i2 
              ENDDO ! i1
              !IF(iw==iw2) print*,'iw dist ',iw,dist(63)
              !
              ! Sort the 125 vectors R by increasing value of |r-R|^2
              ind(1) = 0 ! required for hpsort_eps (see the subroutine)
              CALL hpsort_eps_epw( 125, dist, ind, eps8)
              !
              ! Find all the vectors R with the (same) smallest |r-R|^2;
              ! if R=0 is one of them, then the current point r belongs to 
              ! Wignez-Seitz cell => set found to true
              !
              found = .false.
              i = 1
              mindist = dist(1)
              DO WHILE ( abs(dist(i)-mindist) < eps8 .and. i < 125 )
                IF (ind(i) == 63) found = .true.
                i = i + 1
              ENDDO
              !
              IF (found) THEN
                nrr_tmp(iw,iw2) = nrr_tmp(iw,iw2) + 1
                ndegen_tmp (nrr_tmp(iw,iw2),iw,iw2) = i - 1
                irvec_tmp (:, nrr_tmp(iw,iw2),iw,iw2) = (/ n1, n2, n3 /) 
              ENDIF
            ENDDO ! n3
          ENDDO ! n2 
        ENDDO ! n1 
        !print*,'iw iw2 ',iw, iw2, nrr_tmp(iw,iw2)
      ENDDO ! iw2
    ENDDO ! iw 
    ! 
    nrr = nrr_tmp(1,1)
    irvec(:,:) = irvec_tmp(:,:,1,1) 
    DO iw = 1, dims
      DO iw2 = 1, dims
        DO ir=1, nrr_tmp(iw, iw2)
          found = .false.
          DO irtot = 1, nrr
            IF ( ALL(irvec_tmp(:, ir, iw, iw2) == irvec(:, irtot)) ) THEN
              found = .true.
            ENDIF
          ENDDO !nrr
          IF(.not. found) THEN
            nrr = nrr + 1
            irvec(:, nrr) = irvec_tmp(:, ir, iw, iw2)
          ENDIF
        ENDDO ! ir 
      ENDDO ! nb 
    ENDDO ! na
    ! 
    ! Creates a pair of atoms depended degeneracy array but with a number of WS 
    ! vectors per pair that is equal to the global set. Populate with zero weights 
    ! the one that are not part of that pair set. 
    ndegen(:,:,:) = 0
    DO iw = 1, dims
      DO iw2 = 1, dims
        DO ir=1, nrr_tmp(iw,iw2)
          DO irtot = 1, nrr
            IF ( ALL(irvec(:, irtot) == irvec_tmp(:, ir, iw, iw2)) ) THEN
              ndegen(irtot, iw, iw2) = ndegen_tmp(ir, iw, iw2)
            ENDIF
          ENDDO
        ENDDO
      ENDDO
    ENDDO 
    ! 
    DO iw = 1, dims
      DO iw2 = 1, dims
        tot = 0.d0
        tot2 = 0.d0
        DO i = 1, nrr
          IF (ndegen(i, iw, iw2) > 0 ) THEN
            tot2 = tot2 + 1.d0 / dble(ndegen(i, iw, iw2))
          ENDIF
        ENDDO
        DO i = 1, nrr_tmp(iw, iw2)
          tot = tot + 1.d0 / dble(ndegen_tmp(i, iw, iw2))
        ENDDO
        !
        !print*,'na, nb, tot tot2 ',na,nb,tot,tot2,nq1,nq2,nq3,dble(nq1*nq2*nq3)
        IF(abs(tot-dble(nc1*nc2*nc3)) > eps8) call errore &
         ('wigner_seitzkq',' weights do not add up to nc1*nc2*nc3',1)
        IF(abs(tot-tot2) > eps8) call errore &
         ('wigner_seitzkq',' weigths of pair of atoms is not equal to global weights',1)
      ENDDO
    ENDDO
    ! 
    IF(nrr > 20*nc1*nc2*nc3) call errore &
      ('wigner_seitzkq','too many WS points, try to increase the bound 20*nc1*nc2*nc3',1)
    !
    wslen(:) = 0.d0
    DO i = 1, nrr
      DO ipol = 1, 3
        DO jpol = 1, 3
          wslen(i) = wslen(i) + dble(irvec(ipol,i))*adot(ipol,jpol)*dble(irvec(jpol,i))
        ENDDO
      ENDDO
      wslen(i) = sqrt(wslen(i))
    ENDDO
    !
    CALL cryst_to_cart(dims,shift(:,:),at,1)
    ! 
    !print*,'at/Wan index 1,1'
    !DO i = 1, nrr
    !  if (ndegen(i,1,1)>0) print*,i, irvec(:,i), ndegen(i,1,1)
    !ENDDO
    !print*,'at/Wan index 1,2'
    !DO i = 1, nrr
    !  if (ndegen(i,1,2)>0) print*, i,irvec(:,i), ndegen(i,1,2)
    !ENDDO
    ! 
    DEALLOCATE(ind) 
    ! 
    !-----------------------------------------------------------------------------
    END SUBROUTINE wigner_seitzkq
    !----------------------------------------------------------------------------- 
    ! 
    !-----------------------------------------------------------------
    SUBROUTINE wigner_seitzg (nc1, nc2, nc3, irvec, ndegen, wslen, nrr, w_centers, tau, dims, nat)
    !-----------------------------------------------------------------
    !!
    !! Calculates a grid of points that fall inside of (and eventually 
    !! on the surface of) the Wigner-Seitz supercell centered on 
    !! each atoms. 
    !! Follows Eq. 66 of PRB 55, 10355 (1997).  
    !! We are part of the WS if $R_b - \tau_{\kappa}$ is inside the 
    !! supercell.  
    !! 
    !-----------------------------------------------------------------
    USE kinds,         ONLY : DP
    USE cell_base,     ONLY : at, bg
    USE constants_epw, ONLY : eps8
    !
    implicit none
    !
    INTEGER, INTENT (in) :: nc1
    !! size of the uniform k mesh
    INTEGER, INTENT (in) :: nc2
    !! size of the uniform k mesh
    INTEGER, INTENT (in) :: nc3
    !! size of the uniform k mesh
    INTEGER, INTENT (in) :: dims
    !! Dims is either nbndsub or 1 depending on use_ws
    INTEGER, INTENT (in) :: nat
    !! Number of atoms
    INTEGER, INTENT (out) :: irvec(3,20*nc1*nc2*nc3)
    !! integer components of the ir-th Wigner-Seitz grid point in the basis of the lattice vectors
    INTEGER, INTENT (out) :: ndegen(20*nc1*nc2*nc3,nat,dims,dims)
    !! Number of degeneracies
    INTEGER, INTENT (out) :: nrr
    !! number of Wigner-Seitz grid points 
    REAL(kind=DP), INTENT (in) :: w_centers(3,dims)
    !! Wannier centers
    REAL(kind=DP), INTENT (in) :: tau(3,nat)
    !! Atomic positions
    REAL(kind=DP), INTENT (out) :: wslen(20*nc1*nc2*nc3)
    !! real-space length, in units of alat
    ! 
    ! work variables
    INTEGER :: n1, n2, n3
    !! Index for the larger supercell
    INTEGER :: i1, i2, i3
    !! Index to compute |r-R| distance
    INTEGER :: na, iw, iw2
    !! Atom index
    INTEGER :: i
    !! Iterative index
    INTEGER :: ir
    !! Iterative index on the pair of atoms
    INTEGER :: irtot
    !! Iterative index on the combined pair of atoms
    INTEGER :: ipol, jpol
    !! Cartesian direction 
    INTEGER, ALLOCATABLE :: ind(:)
    !! Index of sorting
    INTEGER :: nind
    !! The metric tensor
    INTEGER :: nrr_tmp(nat,dims,dims)
    !! Temporary array that contains the max number of WS vectors
    !! for a pair of atoms. 
    INTEGER :: irvec_tmp(3,20*nc1*nc2*nc3,nat,dims,dims)
    !! Temporary WS vectors for each atoms 
    INTEGER :: ndegen_tmp(20*nc1*nc2*nc3,nat,dims,dims)
    !! Temporary WS vectors weigths for each atoms
    ! 
    REAL(kind=DP) :: adot(3,3)
    !! Dot product between lattice vector
    REAL(kind=DP) :: dist(125)
    !! Contains the distance squared |r-R|^2
    REAL(kind=DP) :: mindist
    !! Minimum distance
    REAL(kind=DP) :: tot, tot2
    !! Sum of all the weigths
    REAL(kind=DP) :: ndiff(3)
    !! Distances. Must be now real because of fractional atoms
    ! 
    LOGICAL :: found
    !! True if the vector has been found
    ! 
    nind = 20*nc1*nc2*nc3
    IF (nind < 125) THEN
      nind = 125
    ENDIF
    ALLOCATE (ind(nind))
    !
    DO ipol = 1, 3
      DO jpol = 1, 3
        adot (ipol, jpol) = dot_product ( at(:,ipol), at(:,jpol) )
      ENDDO
    ENDDO
    ! 
    CALL cryst_to_cart(nat,tau(:,:),bg,-1)
    CALL cryst_to_cart(dims,w_centers(:,:),bg,-1)
    !
    ! Loop over grid points r on a unit cell that is 8 times larger than a 
    ! primitive supercell. In the end nrr contains the total number of grids 
    ! points that have been found in the Wigner-Seitz cell
    !
    nrr_tmp(:,:,:) = 0
    DO iw = 1, dims
     DO iw2 = 1, dims 
      DO na = 1, nat
        DO n1 = -2*nc1, 2*nc1
          DO n2 = -2*nc2, 2*nc2
            DO n3 = -2*nc3, 2*nc3
              !
              ! Loop over the 5^3 = 125 points R. R=0 corresponds to i1=i2=i3=2, or icnt=63
              !
              i = 0
              dist(:) = 0.d0
              DO i1 = -2, 2
                DO i2 = -2, 2
                  DO i3 = -2, 2
                    i = i + 1
                    !
                    ! Calculate distance squared |r-R|^2 
                    !
                    !ndiff(1) = n1 - i1*nc1 + ( tau(1,na) + w_centers(1,iw2) + w_centers(1,iw) ) / 3.d0
                    !ndiff(2) = n2 - i2*nc2 + ( tau(2,na) + w_centers(2,iw2) + w_centers(2,iw) ) / 3.d0
                    !ndiff(3) = n3 - i3*nc3 + ( tau(3,na) + w_centers(3,iw2) + w_centers(3,iw) ) / 3.d0
                    ndiff(1) = n1 - i1*nc1 + tau(1,na) + ( w_centers(1,iw2) + w_centers(1,iw) ) / 2.d0
                    ndiff(2) = n2 - i2*nc2 + tau(2,na) + ( w_centers(2,iw2) + w_centers(2,iw) ) / 2.d0
                    ndiff(3) = n3 - i3*nc3 + tau(3,na) + ( w_centers(3,iw2) + w_centers(3,iw) ) / 2.d0
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
              ind(1) = 0 ! required for hpsort_eps (see the subroutine)
              CALL hpsort_eps_epw( 125, dist, ind, eps8)
              !
              ! Find all the vectors R with the (same) smallest |r-R|^2;
              ! if R=0 is one of them, then the current point r belongs to 
              ! Wignez-Seitz cell => set found to true
              !
              found = .false.
              i = 1
              mindist = dist(1)
              DO WHILE ( abs(dist(i)-mindist) < eps8 .and. i < 125 )
                IF (ind(i) == 63) found = .true.
                i = i + 1
              ENDDO
              !
              IF (found) THEN
                nrr_tmp(na,iw,iw2) = nrr_tmp(na,iw,iw2) + 1
                ndegen_tmp (nrr_tmp(na,iw,iw2), na, iw, iw2) = i - 1
                irvec_tmp (:, nrr_tmp(na,iw,iw2), na, iw, iw2) = (/ n1, n2, n3 /)
              ENDIF
              !
            ENDDO ! n3   
          ENDDO ! n2 
        ENDDO ! n3
      ENDDO ! na 
     ENDDO ! iw2
    ENDDO ! iw
    ! 
    ! Now creates a global set of WS vectors from all the atoms pair. 
    ! Also remove the duplicated ones. 
    nrr = nrr_tmp(1,1,1)
    irvec(:,:) = irvec_tmp(:, :, 1, 1, 1)
    DO iw = 1, dims
     DO iw2 = 1, dims
      DO na = 1, nat
        DO ir=1, nrr_tmp(na,iw,iw2)
          found = .false.
          DO irtot = 1, nrr
            IF ( ALL(irvec_tmp(:, ir, na, iw, iw2) == irvec(:, irtot)) ) THEN
              found = .true.
            ENDIF
          ENDDO !nrr
          IF(.not. found) THEN
            nrr = nrr + 1
            irvec(:, nrr) = irvec_tmp(:, ir, na, iw, iw2)
          ENDIF
        ENDDO ! ir 
      ENDDO ! na
     ENDDO ! iw2
    ENDDO ! iw
    ! 
    ! Creates a pair of atoms-dependent degeneracy array but with a number of WS 
    ! vectors per pair that is equal to the global set. Populate with zero weights 
    ! the one that are not part of that pair set. 
    ndegen(:,:,:,:) = 0
    DO iw = 1, dims
     DO iw2 = 1, dims
      DO na = 1, nat
        DO ir=1, nrr_tmp(na,iw,iw2)
          DO irtot = 1, nrr
            IF ( ALL(irvec(:, irtot) == irvec_tmp(:, ir, na, iw, iw2)) ) THEN
              ndegen(irtot, na, iw, iw2) = ndegen_tmp(ir, na, iw, iw2)
            ENDIF
          ENDDO
        ENDDO
      ENDDO
     ENDDO
    ENDDO
    ! 
    DO iw = 1, dims
     DO iw2 = 1, dims
      DO na = 1, nat
        tot = 0.d0
        tot2 = 0.d0
        DO i = 1, nrr
          IF (ndegen(i, na, iw, iw2) > 0 ) THEN
            tot2 = tot2 + 1.d0 / dble(ndegen(i, na, iw, iw2))
          ENDIF
        ENDDO
        DO i = 1, nrr_tmp(na, iw, iw2)
          tot = tot + 1.d0 / dble(ndegen_tmp(i, na, iw, iw2))
        ENDDO
        !
        IF(abs(tot-dble(nc1*nc2*nc3)) > eps8) CALL errore &
         ('wigner_seitzg',' weights do not add up to nq1*nq2*nq3',1)
        !WRITE(stdout,*) "tot ", tot, " tot2 ", tot2
        IF(abs(tot-tot2) > eps8) CALL errore &
         ('wigner_seitzg',' weigths of pair of atoms is not equal to global weights',1)
      ENDDO
     ENDDO
    ENDDO
    !
    IF(nrr > 20*nc1*nc2*nc3) CALL errore &
      ('wigner_seitzg','too many WS points, try to increase the bound 20*nc1*nc2*nc3',1)
    !
    ! Now we compute the WS distances
    wslen(:) = 0.d0
    DO i = 1, nrr
      DO ipol = 1, 3
       DO jpol = 1, 3
          wslen(i) = wslen(i) + dble(irvec(ipol,i))*adot(ipol,jpol)*dble(irvec(jpol,i))
       ENDDO
      ENDDO
      wslen(i) = sqrt(wslen(i))
    ENDDO
    !
    !print*,'index band,band,at = 1,1,1'
    !DO i = 1, nrr
    !  if (ndegen(i,1,1,1)>0) print*,i,irvec(:,i), ndegen(i,1,1,1)
    !ENDDO
    !print*,'index band,band,at = 1,2,1'
    !DO i = 1, nrr
    !  if (ndegen(i,1,2,1)>0) print*, i,irvec(:,i), ndegen(i,1,2,1)
    !ENDDO
    !print*,'index band,band,at = 1,2,2'
    !DO i = 1, nrr
    !  if (ndegen(i,1,2,2)>0) print*,i, irvec(:,i), ndegen(i,1,2,2)
    !ENDDO
    !
    CALL cryst_to_cart(nat,tau(:,:),at,1)
    CALL cryst_to_cart(dims,w_centers(:,:),at,1)
    !
    DEALLOCATE(ind)
    ! 
    ! -----------------------------------------------------------------------------------------
    END SUBROUTINE wigner_seitzg
    ! -----------------------------------------------------------------------------------------
    !
  END MODULE wigner
