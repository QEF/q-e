  !                                                                            
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  !--------------------------------------------------------------------------
  SUBROUTINE ephbloch2wanp ( nbnd, nmodes, xk, nq, irvec, wslen, &
    nrk, nrr, epmatwe )
  !--------------------------------------------------------------------------
  !
  !  From the EP Matrix in Electron Bloch representation (coarse mesh), 
  !  find the corresponding matrix in Phonon Wannier representation 
  !
  !--------------------------------------------------------------------------------
  !
  USE kinds,         ONLY : DP
  USE pwcom,         ONLY : at, bg, celldm
  USE elph2,         ONLY : epmatwp
  USE constants_epw, ONLY : bohr2ang, twopi, ci, czero
#ifdef __PARA
  USE io_global,     ONLY : ionode_id
  USE mp,            ONLY : mp_barrier
  USE mp_world,      ONLY : mpime
#endif
  implicit none
  !
  !  input variables - note irvec is dimensioned with nrr_k 
  !                    (which is assumed to be larger than nrr_q)
  !
  integer :: nbnd, nrk, nmodes, nq, nrr, irvec (3, nrk)
  ! number of electronic bands
  ! number of electronic WS points
  ! number of branches
  ! number of qpoints
  ! number of WS points and coordinates
  complex(kind=DP) :: epmatwe (nbnd, nbnd, nrk, nmodes, nq)
  ! EP matrix in electron-wannier representation and phonon bloch representation
  !   (Cartesian coordinates)
  real(kind=DP) :: xk (3, nq), wslen (nrr) 
  ! kpoint coordinates (cartesian in units of 2piba)
  ! WS vectors length (alat units)
  !
  !  output variables
  !
  ! EP matrix in electron-wannier representation and phonon-wannier
  ! representation
  !
  ! work variables 
  !
  integer :: ik, ir, ire
  real(kind=DP) :: rdotk, tmp, rvec1(3), rvec2(3), len1, len2
  complex(kind=DP) :: cfac
  !
  !----------------------------------------------------------
  !  Fourier transform to go into Wannier basis
  !----------------------------------------------------------
  !
  !  D (R) = (1/nk) sum_k e^{-ikR} D (k)
  !
  ! bring xk in crystal coordinates
  !
  CALL cryst_to_cart (nq, xk, at, -1)
  !
  epmatwp = czero
  ! 
  DO ir = 1, nrr
    !
    !
    DO ik = 1, nq
       !
       rdotk = twopi * dot_product( xk ( :, ik), dble(irvec( :, ir) ))
       cfac = exp( -ci*rdotk ) / dble(nq)
       epmatwp(:,:,:,:,ir) = epmatwp(:,:,:,:,ir) + cfac * epmatwe(:,:,:,:,ik)
       !
    ENDDO
    !
    !  check spatial decay of e-p matrix elements in wannier basis - electrons
    !  + phonons
    !
    !  we plot: R_e, R_p, max_{m,n,nu} |g(m,n,nu;R_e,R_p)|
    !
#ifdef __PARA
    IF (mpime.eq.ionode_id) THEN
#endif
      IF (ir.eq.1) open(unit=303,file='decay.epmat_wanep',status='unknown')
      DO ire = 1, nrk
        !
        rvec1 = dble(irvec(1,ire))*at(:,1) + &
                dble(irvec(2,ire))*at(:,2) + &
                dble(irvec(3,ire))*at(:,3)
        rvec2 = dble(irvec(1,ir))*at(:,1) + &
                dble(irvec(2,ir))*at(:,2) + &
                dble(irvec(3,ir))*at(:,3)
        len1 = sqrt(rvec1(1)**2.d0+rvec1(2)**2.d0+rvec1(3)**2.d0)
        len2 = sqrt(rvec2(1)**2.d0+rvec2(2)**2.d0+rvec2(3)**2.d0)
        tmp =  maxval ( abs( epmatwp (:, :, ire, :, ir) ) )
        !
        ! rvec1 : electron-electron0 distance
        ! rvec2 : phonon - electron0 distance
        !
        WRITE(303, '(5f15.10)') len1 * celldm (1) * bohr2ang, &
                                len2 * celldm (1) * bohr2ang, tmp
      ENDDO
      IF (ir.eq.nrr) close(303)
#ifdef __PARA
    ENDIF
#endif
    !
  ENDDO
  !
  ! bring xk back into cart coord
  !
  CALL cryst_to_cart (nq, xk, bg, 1)
  !
  END SUBROUTINE ephbloch2wanp
  !
  ! -----------------------------------------------------------
  !--------------------------------------------------------------------------
  SUBROUTINE ephbloch2wanp_mem ( nbnd, nmodes, xk, nq, irvec, wslen, &
    nrk, nrr, epmatwe )
  !--------------------------------------------------------------------------
  !
  !  From the EP Matrix in Electron Bloch representation (coarse mesh), 
  !  find the corresponding matrix in Phonon Wannier representation 
  !
  !--------------------------------------------------------------------------------
  !
  USE kinds,         ONLY : DP
  USE pwcom,         ONLY : at, bg, celldm
  USE constants_epw, ONLY : bohr2ang, twopi, ci, czero
  USE io_epw,        ONLY : iunepmatwe, iunepmatwp
#ifdef __PARA
  USE io_global,     ONLY : ionode_id
  USE mp,            ONLY : mp_barrier
  USE mp_world,      ONLY : mpime
#endif
  implicit none
  !
  !  input variables - note irvec is dimensioned with nrr_k 
  !                    (which is assumed to be larger than nrr_q)
  !
  integer :: nbnd, nrk, nmodes, nq, nrr, irvec (3, nrk)
  ! number of electronic bands
  ! number of electronic WS points
  ! number of branches
  ! number of qpoints
  ! number of WS points and coordinates
  complex(kind=DP) :: epmatwe (nbnd, nbnd, nrk, nmodes)
  ! EP matrix in electron-wannier representation and phonon bloch representation
  !   (Cartesian coordinates)
  real(kind=DP) :: xk (3, nq), wslen (nrr)
  ! kpoint coordinates (cartesian in units of 2piba)
  ! WS vectors length (alat units)
  !
  !  output variables
  !
  ! EP matrix in electron-wannier representation and phonon-wannier
  ! representation
  !
  !
  ! work variables 
  !
  integer :: ik, ir, ire
  real(kind=DP) :: rdotk, tmp, rvec1(3), rvec2(3), len1, len2
  complex(kind=DP) :: cfac
  COMPLEX(KIND=DP), ALLOCATABLE :: epmatwp_mem (:,:,:,:) !  e-p matrix  in wannier basis 
  !
  !
  ALLOCATE (epmatwp_mem ( nbnd, nbnd, nrk, nmodes))
  ! 
  !----------------------------------------------------------
  !  Fourier transform to go into Wannier basis
  !----------------------------------------------------------
  !
  !  D (R) = (1/nk) sum_k e^{-ikR} D (k)
  !
  ! bring xk in crystal coordinates
  !
  CALL cryst_to_cart (nq, xk, at, -1)
  !
  DO ir = 1, nrr
    !
    epmatwp_mem = czero
    ! 
    DO ik = 1, nq
       !
       ! direct read of epmatwe for this iq 
       CALL rwepmatw ( epmatwe, nbnd, nrk, nmodes, ik, iunepmatwe, -1)
       !
       rdotk = twopi * dot_product( xk ( :, ik), dble(irvec( :, ir) ))
       cfac = exp( -ci*rdotk ) / dble(nq)
       epmatwp_mem = epmatwp_mem + cfac * epmatwe
       !
    ENDDO
    !
    ! direct write of epmatwp_mem for this ir 
    CALL rwepmatw (epmatwp_mem, nbnd, nrk, nmodes, ir, iunepmatwp, +1)
    !  check spatial decay of e-p matrix elements in wannier basis - electrons
    !  + phonons
    !
    !  we plot: R_e, R_p, max_{m,n,nu} |g(m,n,nu;R_e,R_p)|
    !
#ifdef __PARA
    IF (mpime.eq.ionode_id) THEN
#endif
      IF (ir.eq.1) open(unit=303,file='decay.epmat_wanep',status='unknown')
      DO ire = 1, nrk
        !
        rvec1 = dble(irvec(1,ire))*at(:,1) + &
                dble(irvec(2,ire))*at(:,2) + &
                dble(irvec(3,ire))*at(:,3)
        rvec2 = dble(irvec(1,ir))*at(:,1) + &
                dble(irvec(2,ir))*at(:,2) + &
                dble(irvec(3,ir))*at(:,3)
        len1 = sqrt(rvec1(1)**2.d0+rvec1(2)**2.d0+rvec1(3)**2.d0)
        len2 = sqrt(rvec2(1)**2.d0+rvec2(2)**2.d0+rvec2(3)**2.d0)
        tmp =  maxval ( abs( epmatwp_mem (:, :, ire, :) ) )
        !
        ! rvec1 : electron-electron0 distance
        ! rvec2 : phonon - electron0 distance
        !
        WRITE(303, '(5f15.10)') len1 * celldm (1) * bohr2ang, &
                                len2 * celldm (1) * bohr2ang, tmp
      ENDDO
      IF (ir.eq.nrr) close(303)
#ifdef __PARA
    ENDIF
#endif
    !
  ENDDO
  !
  ! bring xk back into cart coord
  !
  CALL cryst_to_cart (nq, xk, bg, 1)
  !
  IF ( ALLOCATED (epmatwp_mem) ) DEALLOCATE (epmatwp_mem)
  !
  END SUBROUTINE ephbloch2wanp_mem





