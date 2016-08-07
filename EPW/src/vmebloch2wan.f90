  !                                                                            
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino 
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  !
  !--------------------------------------------------------------------------
  subroutine vmebloch2wan ( nbnd, nbndsub, nks, nkbl, xk, cu, &
     nrr, irvec, wslen )
  !--------------------------------------------------------------------------
  !!
  !!  Calculate the velocity matrix elements in the Wannier basis
  !!  at no point do we actually have the coarse mesh v-ME. 
  !!
  !--------------------------------------------------------------------------
  !
  USE kinds,     ONLY : DP
  use pwcom,     ONLY : at, bg, celldm, nkstot
  use elph2,     ONLY : cvmew
  USE constants_epw, ONLY : bohr2ang, twopi, zero, ci, czero
  use io_files,  ONLY : prefix
  USE io_global, ONLY : ionode_id
  USE mp_global, ONLY : inter_pool_comm
  USE mp,        ONLY : mp_barrier,mp_sum
  USE mp_world,  ONLY : mpime
  implicit none
  !
  !  input variables
  !
  INTEGER, INTENT (in) :: nbnd
  !! Number of bands
  integer :: nbndsub, nks, nkbl, nrr, irvec (3, nrr), &
       ipool, nkb, nkb_abs, ipol, nnb, ib
  ! number of bands 
  ! number of bands in the optimal subspace 
  ! number of kpoints
  ! number of kpoint blocks, in the pool
  ! number of kpoint blocks, total 
  ! number of WS points and coordinates
  real(kind=DP) :: xk (3, nks), wslen (nrr) , b_tmp(3)
  ! hamiltonian eigenvalues, coarse mesh
  ! kpoint coordinates (cartesian in units of 2piba)
  ! WS vectors length (alat units)
  complex(kind=DP) :: cu (nbnd, nbndsub, nks)
  complex(kind=DP) :: cu_big (nbnd, nbndsub, nkstot)
  complex(kind=DP) :: cfac
  ! rotation matrix from wannier code
  !
  !  output variables
  !
  ! work variables 
  !
  integer :: ik, ibnd, ir, ios, ikstart, ikstop, nkk_abs, nkk
  real(kind=DP) :: rdotk, tmp, zero_vect(3)
  complex(kind=DP) :: Apos(3, nbndsub, nbndsub, nks)
  complex(kind=DP), allocatable :: M_mn(:,:,:,:)
  real(kind=DP), allocatable :: bvec(:,:,:), wb(:)
  character (len=256) :: tempfile
  ! step -1:
  ! setup rotation matrix - we need access to all for the k+b
  cu_big = czero
  CALL ckbounds(ikstart, ikstop)
  cu_big(:,:,ikstart:ikstop) = cu(:,:,:)
  CALL mp_sum(cu_big,inter_pool_comm)
  !
  !  Step 0:
  !  Read in wb, b-vectors
  !
  zero_vect = 0.d0
  tempfile='jesse.vmedat'
  open(1, file=tempfile, action='read', iostat=ios)
  IF (ios /= 0) then
     !
     ! end up leaving zeros for everything.  In this case
     ! obviously the velocities will be meaningless
     ! This should allow the program to run, however
     !
     nnb = 1
     allocate (  M_mn(nbnd, nbnd, nnb, nkstot),  &
          bvec(3,nnb,nkstot),             &
          wb(nnb) )
     bvec = zero
     wb   = zero
  ELSE
     read(1,*) nnb
     allocate (  M_mn(nbnd, nbnd, nnb, nkstot),  &
                 bvec(3,nnb,nkstot),             &
                 wb(nnb) )
     DO ik = 1, nkstot
        DO ib = 1, nnb
           read(1,'(4f20.10)') bvec(:,ib,ik), wb(ib)
        ENDDO
     ENDDO
     close(1)
  ENDIF
  !
  bvec = bvec * bohr2ang
  wb = wb / bohr2ang**2
  !  read Mmn for velocity calculation
  tempfile=trim(prefix)//'.mmn'
  open(1, file=tempfile, form='unformatted', action='read', iostat=ios)
  IF (ios /= 0) then
     ! if it doesn't exist, then we just set the mmn to zero.  I'll have to 
     ! clean this up
     CALL errore ('vmebloch2wan','error opening' // tempfile, 0)
     M_mn = czero
  ELSE
     !
     read(1) M_mn
     close(1)
     !
  ENDIF
  !
  !  Step 0.1
  ! Calculate (<u_kn^W | u_(k+b)m^W> - delta_mn)
  !
  DO ik = 1, nks
     DO ib = 1, nnb
        !
        CALL ktokpmq ( xk(:,ik),zero_vect, +1,ipool,nkk,nkk_abs)
        b_tmp(:) = celldm(1) /(twopi) * bvec(:,ib,nkk_abs)
        CALL ktokpmq ( xk(:,ik),b_tmp(:), +1,ipool,nkb,nkb_abs)
        !
        M_mn (:,:, ib, ik ) = matmul (  conjg(transpose(cu(:,:,ik))), M_mn(:,:,ib,nkk_abs) )
        M_mn (:,:, ib, ik ) = matmul ( M_mn(:,:,ib,ik), cu_big(:,:,nkb_abs) )
        !
        DO ibnd = 1, nbndsub
           M_mn(ibnd,ibnd, ib, ik) = M_mn(ibnd,ibnd, ib, ik)  - 1.d0
        ENDDO
        !
     ENDDO
  ENDDO
  !
  ! Calculate A_mn(k)^(W) [Eqn. 44 of PRB 74 195118 (2006)]
  !
  Apos = czero
  DO ik = 1, nks
     CALL ktokpmq ( xk(:,ik),zero_vect, +1,ipool,nkk,nkk_abs)
     !
     DO ib =1, nnb
        DO ipol = 1, 3
           Apos( ipol, :, :, ik) = Apos (ipol, :,:,ik) + &
                ci * wb(ib) * bvec(ipol,ib, nkk_abs) * M_mn(1:nbndsub,1:nbndsub,ib, ik)
        ENDDO
     ENDDO
  ENDDO
  !
  !----------------------------------------------------------
  !  Fourier transform to go into Wannier basis
  !----------------------------------------------------------
  !
  !
  ! bring xk in crystal coordinates
  !
  CALL cryst_to_cart (nks, xk, at, -1)
  !
  cvmew ( :, :, :, :) = czero
  ! 
  DO ir = 1, nrr
    !
    DO ik = 1, nks
       !
       rdotk = twopi * dot_product( xk ( :, ik), dble(irvec( :, ir) ))
       cfac = exp( -ci*rdotk ) / dble(nkbl)
       cvmew ( :, :, :, ir ) = cvmew ( :, :, :, ir ) + cfac * Apos ( :, :, :, ik )
       !
    ENDDO
    !
  ENDDO
  CALL mp_sum(cvmew,inter_pool_comm) 
  !
  ! bring xk back into cart coord
  !
  CALL cryst_to_cart (nks, xk, bg, 1)
  !
  !
  !  check spatial decay of velocity matrix elements in Wannier basis
  !  the unit in r-space is angstrom
  !
  IF (mpime.eq.ionode_id) then
    open(unit=300,file='decay.v')
    WRITE(300, '(/3x,a/)') '#Spatial decay of Velocity matrix element in Wannier basis'
    DO ir = 1, nrr
      !
      tmp =  maxval ( abs( cvmew (:,:,:,ir)) )
      WRITE(300, *) wslen(ir) * celldm (1) * bohr2ang, tmp
      !
    ENDDO
    close(300)
  ENDIF
  CALL mp_barrier(inter_pool_comm)
  !
  WRITE(6,'(/5x,a)') 'Velocity matrix elements calculated'
  WRITE(6,*)
  !
  end subroutine vmebloch2wan
