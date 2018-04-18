  !
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino  
  ! 
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .
  !
  !----------------------------------------------------------------------
  MODULE bloch2wan
  !----------------------------------------------------------------------
  !!  
  !! Contains all the routines that transforms quantities from Bloch space 
  !! to Wannier space. 
  !! 
  !!  
  IMPLICIT NONE
  ! 
  CONTAINS
    !
    !--------------------------------------------------------------------------
    SUBROUTINE hambloch2wan ( nbnd, nbndsub, nks, nkstot, et, xk, cu, &
       lwin, nrr, irvec, wslen, chw )
    !--------------------------------------------------------------------------
    !!
    !!  From the Hamiltonian in Bloch representationi (coarse mesh), 
    !!  find the corresponding Hamiltonian in Wannier representation 
    !!
    !--------------------------------------------------------------------------
    !
    USE kinds,     ONLY : DP
    USE cell_base, ONLY : at, bg, celldm
    USE constants_epw, ONLY : bohr2ang, twopi, ci, czero
    USE io_global, ONLY : ionode_id
    USE mp_global, ONLY : inter_pool_comm
    USE mp,        ONLY : mp_barrier,mp_sum
    USE mp_world,  ONLY : mpime
    !
    implicit none
    !
    INTEGER, INTENT (in) :: nbnd
    !! number of bands 
    INTEGER, INTENT (in) :: nks
    !! number of kpoints
    ! 
    LOGICAL, INTENT (in) :: lwin( nbnd, nks )
    !! identify bands within outer energy window (for disentanglement)
    ! 
    INTEGER, INTENT (in) :: nbndsub
    !! number of bands in the optimal subspace 
    INTEGER, INTENT (in) :: nkstot
    !! number of kpoint blocks, in the pool
    INTEGER, INTENT (in) ::  nrr
    !! number of kpoint blocks, total 
    INTEGER, INTENT (in) :: irvec (3, nrr)
    !! number of WS points and coordinates
    !
    REAL(kind=DP), INTENT (in) :: et (nbnd, nks)
    !! hamiltonian eigenvalues, coarse mesh
    REAL(kind=DP), INTENT (in) :: xk (3, nks)
    !! kpoint coordinates (cartesian in units of 2piba)
    REAL(kind=DP), INTENT (in) :: wslen (nrr)
    !! WS vectors length (alat units)
    !
    COMPLEX(kind=DP), INTENT (in) :: cu (nbnd, nbndsub, nks)
    !! rotation matrix from wannier code
    COMPLEX(kind=DP), INTENT(OUT) :: chw ( nbndsub, nbndsub, nrr)
    !! Hamiltonian in smooth Bloch basis, coarse mesh 
    !
    ! work variables 
    !
    INTEGER :: ik
    !! Counter of k-point index
    INTEGER :: ibnd
    !! Counter on band index
    INTEGER :: jbnd
    !! Counter on band index
    INTEGER :: ir
    !! Counter on real-space index
    INTEGER :: mbnd
    !! Counter on band index 
    INTEGER :: j
    !! Counter on band index
    REAL(kind=DP) :: rdotk
    !! $$\mathbf{r}\cdot\mathbf{k}
    REAL(kind=DP) :: tmp
    !! Maximum value of the real space Hamiltonian
    REAL(kind=dp)    :: et_opt(nbnd,nks)
    !! hamiltonian eigenvalues within the outer window in the first ndimwin(ik) entries
    !
    COMPLEX(kind=DP) :: chs(nbndsub, nbndsub, nks) 
    !! Hamiltonian
    COMPLEX(kind=DP) :: cfac
    !! $$e^{-i\mathbf{r}\cdot\mathbf{k}}$$
    COMPLEX(kind=DP) :: ctmp
    !! Temporary variable to store the Hamiltonian
    !
    !----------------------------------------------------------
    !    STEP 1: rotation to optimally smooth Bloch states
    !----------------------------------------------------------
    !
    CALL start_clock ( 'Ham: step 1' )
    !
    et_opt=0.d0
    !
    ! slim down et to contain states within the outer window
    !
    do ik=1,nks
       ibnd=0
       do j=1,nbnd
          if(lwin(j,ik)) then
             ibnd=ibnd+1
             et_opt(ibnd,ik)=et(j,ik)
          end if
       end do
    end do
    !
    !  H~ (k) = U(k)^\dagger * H(k) * U(k)
    !  H~ (k) is chs( nbndsub, nbndsub, ik )
    !
    DO ik = 1, nks
     !
     !
     DO jbnd = 1, nbndsub
      DO ibnd = 1, jbnd
         !
         ctmp = czero
         !
         DO mbnd = 1, nbnd
           ctmp = ctmp + conjg(cu (mbnd,ibnd,ik)) * et_opt (mbnd,ik) * cu (mbnd,jbnd,ik)
         ENDDO
         !
         chs (ibnd , jbnd , ik) = ctmp 
         chs (jbnd , ibnd , ik) = conjg(ctmp)
         !
      ENDDO
     ENDDO
    ENDDO
    !
    CALL stop_clock ( 'Ham: step 1' )
    !
    !----------------------------------------------------------
    !  STEP 2: Fourier transform to go into Wannier basis
    !----------------------------------------------------------
    !
    !  H (R) = (1/nk) sum_k e^{-ikR} H~ (k)
    !  chw (nbndsub, nbndsub, ir) is H (R)
    !
    CALL start_clock ( 'Ham: step 2' )
    !
    ! bring xk in crystal coordinates
    !
    CALL cryst_to_cart (nks, xk, at, -1)
    !
    chw ( :, :, :) = czero 
    ! 
    DO ir = 1, nrr
      !
      DO ik = 1, nks
         !
         rdotk = twopi * dot_product( xk ( :, ik), dble(irvec( :, ir) ))
         cfac = exp( -ci*rdotk ) / dble(nkstot)
         chw ( :, :, ir ) = chw ( :, :, ir ) + cfac * chs ( :, :, ik )
         !
      ENDDO
      !
    ENDDO
    CALL mp_sum(chw,inter_pool_comm) 
    !
    ! bring xk back into cart coord
    !
    CALL cryst_to_cart (nks, xk, bg, 1)
    !
      !
      !  check spatial decay of Hamiltonian in Wannier basis
      !  the unit in r-space is angstrom
      !
      IF (mpime.eq.ionode_id) THEN
         open(unit=300,file='decay.H')
         WRITE(300, '(/3x,a/)') '#Spatial decay of Hamiltonian in Wannier basis'
        DO ir = 1, nrr
          !
          tmp =  maxval ( abs( chw (:,:,ir)) )
          WRITE(300, *) wslen(ir) * celldm (1) * bohr2ang, tmp
          !
        ENDDO
        close(300)
      ENDIF
      CALL mp_barrier(inter_pool_comm)
    !
    CALL stop_clock ( 'Ham: step 2' )
    !
    END SUBROUTINE hambloch2wan 
    !
    !--------------------------------------------------------------------------
    SUBROUTINE dmebloch2wan ( nbnd, nbndsub, nks, nkbl, dmec, xk, cu, &
       nrr, irvec, wslen, lwin )
    !--------------------------------------------------------------------------
    !!
    !!  From the Dipole in Bloch representationi (coarse mesh), 
    !!  find the corresponding Dipole in Wannier representation 
    !!
    !
    USE kinds,         ONLY : DP
    USE cell_base,     ONLY : at, bg, celldm
    USE elph2,         ONLY : cdmew
    USE constants_epw, ONLY : bohr2ang, twopi, ci, czero
    USE io_global,     ONLY : ionode_id
    USE mp_global,     ONLY : inter_pool_comm
    USE mp_world,      ONLY : mpime
    USE mp,            ONLY : mp_barrier,mp_sum
    implicit none
    !
    !  input variables
    !
    INTEGER, INTENT (in) :: nbnd
    !! number of bands
    INTEGER, INTENT (in) :: nbndsub
    !! number of bands in the optimal subspace
    INTEGER, INTENT (in) :: nks
    !! number of kpoints
    INTEGER, INTENT (in) :: nkbl
    !! number of kpoint blocks, in the pool
    INTEGER, INTENT (in) :: nrr
    !! number of WS points 
    INTEGER, INTENT (in) :: irvec (3, nrr) 
    !! Coordinate of Wannier space points
    ! 
    REAL(kind=DP), INTENT (in) :: xk (3, nks)
    !! kpoint coordinates (cartesian in units of 2piba) 
    REAL(kind=DP), INTENT (in) :: wslen (nrr)
    !! WS vectors length (alat units)
    ! 
    COMPLEX(kind=DP), INTENT (in) :: dmec (3,nbnd, nbnd,nks)
    !! Dipole matrix elements on coarse mesh
    COMPLEX(kind=DP), INTENT (in) :: cu (nbnd, nbndsub, nks)
    !! rotation matrix from wannier code
    !
    LOGICAL, INTENT(in) :: lwin(nbnd,nks)
    !! 
    !
    ! Local variables
    INTEGER :: ipol
    !! Counter on polarization
    INTEGER :: ik
    !! Counter on k-point
    INTEGER :: ir
    !! Counter on WS points
    INTEGER :: i
    !! Counter on total band index
    INTEGER :: j
    !! Counter on total band index
    INTEGER :: ibnd
    !! Counter on band that are in the Wannierized window
    INTEGER :: jbnd
    !! Counter on band that are in the Wannierized window
    !
    REAL(kind=DP) :: rdotk
    !! $$ mathbf{r}\cdot\mathbf{k} $$
    REAL(kind=DP) :: tmp
    !! Temporary variables 
    COMPLEX(kind=DP) :: cps(3, nbndsub, nbndsub, nks)
    !! Hamiltonian in smooth Bloch basis, coarse mesh 
    COMPLEX(kind=DP) :: cfac
    !! $$ e^{-i\mathbf{r}\cdot\mathbf{k}} $$
    COMPLEX(kind=DP) :: dmec_utmp(nbnd,nbndsub)
    !! dmec after multiplication with the Wannier rotation matrix cu.
    COMPLEX(kind=DP) :: dmec_opt(3, nbnd,nbnd, nks)
    !! dmec computed in pmn, rescaled down if skipping lwin.
    !
    CALL start_clock ( 'Dipole: step 1' )
    !
    !--------------------------------------------------------------
    !    STEP 0: Rescale the optical matrix on the coarse grid down
    !            This is if you skip band during the Wannierization
    !--------------------------------------------------------------
    ! 
    dmec_opt=czero
    DO ik=1,nks
       ibnd=0
       DO i=1,nbnd
          IF(lwin(i,ik)) THEN
             ibnd=ibnd+1
             jbnd=0
             DO j=1,nbnd
                IF(lwin(j,ik)) THEN
                   jbnd=jbnd+1
                   dmec_opt(:,ibnd,jbnd,ik)=dmec(:,i,j,ik)
                END IF
             END DO
          END IF
       END DO
    END DO
    ! 
    !----------------------------------------------------------
    !    STEP 1: rotation to optimally smooth Bloch states
    !----------------------------------------------------------
    !
    !  p~ (k) = U(k)^\dagger * p(k) * U(k)
    !  p~ (k) is cps( ipol, nbndsub, nbndsub, ik )
    !
    DO ik = 1, nks
       DO ipol = 1, 3
          ! copied from ephbloch2wane.  produce equivalent results
          !
          dmec_utmp(:,:) = matmul( dmec_opt(ipol,:,:,ik), cu(:,:,ik) )
          cps (ipol, :,:, ik) = matmul ( conjg(transpose( cu(:,:,ik))), dmec_utmp (:,:) )
          !
      ENDDO
    ENDDO
    !
    !----------------------------------------------------------
    !  STEP 2: Fourier transform to go into Wannier basis
    !----------------------------------------------------------
    !
    !  p (R) = (1/nk) sum_k e^{-ikR} p~ (k)
    !  cdmew (ipol, nbndsub, nbndsub, ir) is p (R)
    !
    CALL start_clock ( 'Dipole: step 2' )
    !
    ! bring xk in crystal coordinates
    !
    CALL cryst_to_cart (nks, xk, at, -1)
    !
    cdmew ( :, :, :, :) = czero 
    ! 
    DO ir = 1, nrr
      !
      DO ik = 1, nks
         !
         rdotk = twopi * dot_product( xk ( :, ik), dble(irvec( :, ir) ))
         cfac = exp( -ci*rdotk ) / dble(nkbl)
         DO ipol = 1, 3
            cdmew ( ipol, :, :, ir ) = cdmew ( ipol, :, :, ir ) + cfac * cps ( ipol, :, :, ik )
         ENDDO
         !
      ENDDO
      !
    ENDDO
    CALL mp_sum(cdmew,inter_pool_comm) 
    !
    ! bring xk back into cart coord
    !
    CALL cryst_to_cart (nks, xk, bg, 1)
    !
    ! Check spatial decay of Dipole in Wannier basis
    ! the unit in r-space is angstrom
    !
    IF (mpime.eq.ionode_id) then
       OPEN(unit=300,file='decay.P')
       WRITE(300, '(/3x,a/)') '#Spatial decay of Dipole in Wannier basis'
       DO ir = 1, nrr
          !
          tmp =  maxval ( abs( cdmew (:, :,:,ir)) )
          WRITE(300, *) wslen(ir) * celldm (1) * bohr2ang, tmp
          !
       ENDDO
       close(300)
    ENDIF
    CALL mp_barrier(inter_pool_comm)
    !
    CALL stop_clock ( 'Dipole: step 2' )
    !
    END SUBROUTINE dmebloch2wan
    ! 
    !------------------------------------------------------------------------
    SUBROUTINE dynbloch2wan ( nmodes, nq, xk, dynq, nrr, irvec, wslen )
    !------------------------------------------------------------------------
    !!
    !!  From the Dynamical Matrix in Bloch representation (coarse mesh), 
    !!  find the corresponding matrix in Wannier representation 
    !!
    !!  NOTA BENE: it seems to be very important that the matrix is kept real
    !!  physically these are truely the interatomic force constants.
    !!  If you use a complex matrix instead, you may get some spurious 
    !!  oscillations when you interpolate the phonon dispersions.
    !!
    !!  Note also that the acoustic sum rule for the q=0 case has been imposed
    !!  already in readmat_shuffle.f90
    !!
    !
    USE kinds,         ONLY : DP
    USE cell_base,     ONLY : at, bg, celldm
    USE ions_base,     ONLY : nat, tau
    USE phcom,         ONLY : nq1, nq2, nq3
    USE elph2,         ONLY : rdw, epsi, zstar
    USE epwcom,        ONLY : lpolar
    USE constants_epw, ONLY : bohr2ang, twopi, ci
    USE io_global,     ONLY : ionode_id
    USE mp_world,      ONLY : mpime
    USE mp,            ONLY : mp_barrier
    USE mp_global,     ONLY : inter_pool_comm
    ! 
    implicit none
    !
    !  input variables
    !
    integer :: nmodes, nq, nrr, irvec (3, nrr)
    ! number of branches
    ! number of qpoints
    ! number of WS points and coordinates
    complex(kind=DP) :: dynq (nmodes, nmodes, nq)
    ! dynamical matrix in bloch representation (Cartesian coordinates)
    real(kind=DP) :: xk (3, nq), wslen (nrr) 
    ! kpoint coordinates (cartesian in units of 2piba)
    ! WS vectors length (alat units)
    !
    !  output variables
    !
    ! work variables 
    !
    integer :: ik, ir
    real(kind=DP) :: rdotk, tmp
    complex(kind=DP) :: cfac
    !
    !  subtract the long-range term from D(q)
    IF (lpolar) THEN
      DO ik = 1, nq
        !xk has to be in cart. coord.
        CALL rgd_blk (nq1,nq2,nq3,nat,dynq(1,1,ik),xk(:,ik),tau,epsi,zstar,-1.d0)
        !
      ENDDO
    ENDIF
    !
    !----------------------------------------------------------
    !  Fourier transform to go into Wannier basis
    !----------------------------------------------------------
    !
    !  D (R) = (1/nk) sum_k e^{-ikR} D (k)
    !  rdw (nmodes, nmodes, ir) is D (R)
    !
    ! bring xk in crystal coordinates
    !
    CALL cryst_to_cart (nq, xk, at, -1)
    !
    rdw ( :, :, :) = 0.d0
    ! 
    DO ir = 1, nrr
      !
      DO ik = 1, nq
         !
         rdotk = twopi * dot_product( xk ( :, ik), dble(irvec( :, ir) ))
         cfac = exp( -ci*rdotk ) / dble(nq)
      !DBSP - real was commented
         rdw ( :, :, ir ) = rdw ( :, :, ir ) +  cfac * dynq ( :, :, ik ) 
        ! rdw ( :, :, ir ) = rdw ( :, :, ir ) + real ( cfac * dynq ( :, :, ik ) ) 
         !                                     ^^^^
         !                                   note this
         !
      ENDDO
      !
    ENDDO
    !
    ! bring xk back into cart coord
    !
    CALL cryst_to_cart (nq, xk, bg, 1)
    !
    !
    !  check spatial decay of dynamical matrix in Wannier basis
    !  the unit in r-space is angstrom, and I am plotting
    !  the matrix for the first mode only
    !
    IF (mpime.eq.ionode_id) THEN
      OPEN(unit=302,file='decay.dynmat')
      WRITE(302, '(/3x,a/)') '#Spatial decay of Dynamical matrix in Wannier basis'
      DO ir = 1, nrr
        !
        tmp =  maxval ( abs( rdw (:,:,ir)) )
        WRITE(302, *) wslen(ir) * celldm (1) * bohr2ang, tmp
        !
      ENDDO
      CLOSE(302)
    ENDIF
    CALL mp_barrier(inter_pool_comm)
    !
    END SUBROUTINE dynbloch2wan
    !-----------------------------------------------------
    !--------------------------------------------------------------------------
    SUBROUTINE vmebloch2wan ( nbnd, nbndsub, nks, nkbl, xk, cu, &
       nrr, irvec, wslen )
    !--------------------------------------------------------------------------
    !!
    !!  Calculate the velocity matrix elements in the Wannier basis
    !!  at no point do we actually have the coarse mesh v-ME. 
    !!
    !--------------------------------------------------------------------------
    !
    USE kinds,     ONLY : DP
    use cell_base, ONLY : at, bg, celldm
    use pwcom,     ONLY : nkstot
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
    END SUBROUTINE vmebloch2wan
    ! 
    !-----------------------------------------------------------------------
    SUBROUTINE ephbloch2wane ( nbnd, nbndsub, nks, nkstot, xk, &
         cu, cuq, epmatk, nrr, irvec, wslen, epmatw)
    !-----------------------------------------------------------------------
    !!
    !!  From the electron-phonon matrix elements in Bloch representation (coarse 
    !!  mesh), find the corresponding matrix elements in Wannier representation
    !!
    !!
    !-----------------------------------------------------------------------
    !
    USE kinds,     ONLY : DP
    USE cell_base, ONLY : at, bg, celldm
    USE constants_epw, ONLY : bohr2ang, twopi, ci, czero, cone
    USE io_epw,    ONLY : iuwane
    USE io_global, ONLY : ionode_id
    USE mp_global, ONLY : inter_pool_comm
    USE mp       , ONLY : mp_sum 
    USE mp_world,  ONLY : mpime
    implicit none
    !
    !  input variables
    !
    INTEGER, INTENT (in) :: nbnd
    !! Number of bands
    INTEGER, INTENT (in) :: nbndsub
    !! Number of bands in the optimal subspace
    INTEGER, INTENT (in) :: nks
    !! Number of kpoints
    INTEGER, INTENT (in) :: nrr
    !! Number of kpoint blocks, in the pool
    INTEGER, INTENT (in) :: nkstot
    !! Number of kpoint blocks, total
    INTEGER, INTENT (in) :: irvec (3, nrr)
    !! Number of WS points and coordinates
    !
    REAL(kind=DP), INTENT (in) :: xk (3, nks)
    !! kpoint coordinates (cartesian in units of 2piba)
    REAL(kind=DP), INTENT (in) :: wslen (nrr)
    !! WS vectors length (alat units)
    COMPLEX(kind=DP), INTENT (in) :: cu (nbnd, nbndsub, nks)
    !! Rotation matrix from wannier code
    COMPLEX(kind=DP), INTENT (in) :: cuq (nbnd, nbndsub, nks)
    !! Rotation matrix from wannier code
    COMPLEX(kind=DP), INTENT (in) :: epmatk ( nbnd, nbnd, nks)
    !! e-p matrix in bloch representation, coarse mesh
    !
    ! output variables 
    !
    COMPLEX(kind=DP), INTENT(out) :: epmatw ( nbndsub, nbndsub, nrr)
    !!  e-p matrix  in wannier basis 
    !
    ! Work variables 
    !
    complex(kind=DP) :: epmats (nbndsub, nbndsub, nks), eptmp(nbndsub, nbnd)
    !  e-p matrix  in smooth Bloch basis, coarse mesh
    !  e-p matrix, temporary
    !
    integer :: ik, ir
    real(kind=DP) :: rdotk, tmp
    complex(kind=DP) :: cfac
    !
    !
    !----------------------------------------------------------
    !  STEP 1: rotation to optimally smooth Bloch states
    !----------------------------------------------------------
    !
    !  g~ = U_k+q^\dagger g U_k
    !
    !  g   is epmatk (ibnd, jbnd, ik)
    !  g~  is epmats (ibnd, jbnd, ik)
    !
    CALL start_clock ( 'ep: step 1' )
    !
    !
    DO ik = 1, nks
       !
       ! the two zgemm calls perform the following ops:
       ! epmats  = [ cu(ikq)^\dagger * epmatk ] * cu(ikk)
       ! [here we have a size-reduction from nbnd*nbnd to nbndsub*nbndsub] 
       !
       CALL zgemm ('c', 'n', nbndsub, nbnd, nbnd, cone, cuq(:,:,ik),  &
            nbnd, epmatk(:,:,ik), nbnd, czero, eptmp, nbndsub)
       CALL zgemm ('n', 'n', nbndsub, nbndsub, nbnd, cone, eptmp,     &
            nbndsub, cu(:,:,ik), nbnd, czero, epmats(:,:,ik), nbndsub)
       !
    ENDDO
    !
    CALL stop_clock ( 'ep: step 1' )
    !
    !----------------------------------------------------------------------
    !  STEP 2: Fourier transform to obtain matrix elements in wannier basis
    !----------------------------------------------------------------------
    !
    !  g (R) = (1/nkc) sum_k e^{-ikR} g~(k)
    !
    !  epmatw (nbndsub,nbndsub,ir) is g(R)
    !
    CALL start_clock ( 'ep: step 2' )
    !
    epmatw (:, :, :) = czero
    !
    ! bring xk in crystal coordinates
    !
    CALL cryst_to_cart (nks, xk, at, -1)
    !
    DO ir = 1, nrr
       !
       DO ik = 1, nks
         !
         !
         rdotk = twopi * dot_product( xk ( :, ik), dble(irvec( :, ir) ))
         cfac = exp( -ci*rdotk ) / dble(nkstot)
         epmatw ( :, :, ir) = epmatw ( :, :, ir) + cfac * epmats ( :, :, ik)
         !
       ENDDO
       !
    ENDDO
    !
    CALL mp_sum(epmatw,inter_pool_comm)  
    !
    ! bring xk back into cart coord
    !
    CALL cryst_to_cart (nks, xk, bg, 1)
    !
    !
    !  Sheck spatial decay of matrix elements in Wannier basis
    !  the unit in r-space is angstrom, and I am plotting 
    !  the matrix for the first mode only
    !
    IF (mpime.eq.ionode_id) THEN
      OPEN (unit=iuwane,file='decay.epwane')
      WRITE(iuwane, '(a)') '# Spatial decay of e-p matrix elements in Wannier basis'
      DO ir = 1, nrr
        ! 
        tmp =  maxval ( abs(epmatw(:,:,ir)) ) 
        WRITE(iuwane, *) wslen(ir) * celldm (1) * bohr2ang, tmp
        !
      ENDDO
      CLOSE(iuwane)
    ENDIF
    !
    CALL stop_clock ( 'ep: step 2' )
    !
    END SUBROUTINE ephbloch2wane  
    ! 
    !--------------------------------------------------------------------------
    SUBROUTINE ephbloch2wanp ( nbnd, nmodes, xk, nq, irvec, &
      nrk, nrr, epmatwe )
    !--------------------------------------------------------------------------
    !!
    !!  From the EP Matrix in Electron Bloch representation (coarse mesh), 
    !!  find the corresponding matrix in Phonon Wannier representation 
    !!
    !--------------------------------------------------------------------------------
    !
    USE kinds,         ONLY : DP
    USE cell_base,     ONLY : at, bg, celldm
    USE elph2,         ONLY : epmatwp
    USE constants_epw, ONLY : bohr2ang, twopi, ci, czero
    USE io_epw,        ONLY : iuwanep
    USE io_global,     ONLY : ionode_id
    USE mp,            ONLY : mp_barrier
    USE mp_world,      ONLY : mpime
    implicit none
    !
    !  Input variables - note irvec is dimensioned with nrr_k 
    !                    (which is assumed to be larger than nrr_q)
    !
    INTEGER, INTENT(in) :: nbnd
    !! Number of electronic bands
    INTEGER, INTENT(in) :: nrk
    !! number of electronic WS points
    INTEGER, INTENT(in) :: nmodes
    !! number of branches
    INTEGER, INTENT(in) :: nq 
    !! number of qpoints
    INTEGER, INTENT(in) :: nrr
    !! number of WS points and coordinates
    INTEGER, INTENT(in) :: irvec (3, nrk)
    !! Real space vector (irvec is dimensioned with nrr_k)
    !
    REAL(kind=DP), INTENT(in) :: xk (3, nq)
    !! Kpoint coordinates (cartesian in units of 2piba) 
    ! 
    COMPLEX(kind=DP), INTENT(in) :: epmatwe (nbnd, nbnd, nrk, nmodes, nq)
    !! EP matrix in electron-wannier representation and phonon bloch representation
    !!   (Cartesian coordinates)
    !
    !  Output variables
    !
    ! EP matrix in electron-wannier representation and phonon-wannier
    ! representation
    !
    ! Work variables 
    !
    INTEGER :: ik, ir, ire
    REAL(kind=DP) :: rdotk, tmp, rvec1(3), rvec2(3), len1, len2
    COMPLEX(kind=DP) :: cfac
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
      IF (mpime.eq.ionode_id) THEN
        IF (ir.eq.1) open(unit=iuwanep,file='decay.epmat_wanep',status='unknown')
        IF (ir.eq.1) WRITE(iuwanep, '(a)') '#  R_e,    R_p, max_{m,n,nu} |g(m,n,nu;R_e,R_p)| '
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
          WRITE(iuwanep, '(5f15.10)') len1 * celldm (1) * bohr2ang, &
                                  len2 * celldm (1) * bohr2ang, tmp
        ENDDO
        IF (ir.eq.nrr) close(iuwanep)
      ENDIF
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
    SUBROUTINE ephbloch2wanp_mem ( nbnd, nmodes, xk, nq, irvec, &
      nrk, nrr, epmatwe )
    !--------------------------------------------------------------------------
    !
    !  From the EP Matrix in Electron Bloch representation (coarse mesh), 
    !  find the corresponding matrix in Phonon Wannier representation 
    !
    !--------------------------------------------------------------------------------
    !
    USE kinds,         ONLY : DP
    USE cell_base,     ONLY : at, bg, celldm
    USE constants_epw, ONLY : bohr2ang, twopi, ci, czero
    USE io_epw,        ONLY : iunepmatwe, iunepmatwp, iuwanep
    USE io_global,     ONLY : ionode_id
    USE mp,            ONLY : mp_barrier
    USE mp_world,      ONLY : mpime
    implicit none
    !
    !  input variables - note irvec is dimensioned with nrr_k 
    !                    (which is assumed to be larger than nrr_q)
    !
    INTEGER, INTENT(in) :: nbnd
    !! Number of electronic bands
    INTEGER, INTENT(in) :: nrk
    !! number of electronic WS points
    INTEGER, INTENT(in) :: nmodes
    !! number of branches
    INTEGER, INTENT(in) :: nq
    !! number of qpoints
    INTEGER, INTENT(in) :: nrr
    !! number of WS points and coordinates
    INTEGER, INTENT(in) :: irvec (3, nrk)
    !! Real space vector (irvec is dimensioned with nrr_k)
    !
    REAL(kind=DP), INTENT(in) :: xk (3, nq)
    !! Kpoint coordinates (cartesian in units of 2piba) 
    ! 
    COMPLEX(kind=DP), INTENT(in) :: epmatwe (nbnd, nbnd, nrk, nmodes)
    !! EP matrix in electron-wannier representation and phonon bloch representation
    !!   (Cartesian coordinates)
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
      IF (mpime.eq.ionode_id) THEN
        IF (ir.eq.1) open(unit=iuwanep,file='decay.epmat_wanep',status='unknown')
        IF (ir.eq.1) WRITE(iuwanep, '(a)') '#  R_e,    R_p, max_{m,n,nu} |g(m,n,nu;R_e,R_p)| '
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
          WRITE(iuwanep, '(5f15.10)') len1 * celldm (1) * bohr2ang, &
                                  len2 * celldm (1) * bohr2ang, tmp
        ENDDO
        IF (ir.eq.nrr) close(iuwanep)
      ENDIF
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
    ! 
  END MODULE bloch2wan
