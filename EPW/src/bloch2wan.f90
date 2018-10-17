  !
  !k Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino  
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
       lwin, exband, nrr, irvec, wslen, chw )
    !--------------------------------------------------------------------------
    !!
    !!  From the Hamiltonian in Bloch representationi (coarse mesh), 
    !!  find the corresponding Hamiltonian in Wannier representation 
    !!
    !--------------------------------------------------------------------------
    !
    USE kinds,     ONLY : DP
    USE cell_base, ONLY : at, bg, celldm
    USE constants_epw, ONLY : bohr2ang, twopi, ci, czero, zero, ryd2ev
    USE io_global, ONLY : ionode_id
    USE io_epw,    ONLY : iudecayH
    USE mp_global, ONLY : inter_pool_comm
    USE mp,        ONLY : mp_barrier, mp_sum
    USE mp_world,  ONLY : mpime
    !
    IMPLICIT NONE
    !
    LOGICAL, INTENT (in) :: lwin( nbnd, nks )
    !! Bands at k within outer energy window
    LOGICAL, INTENT(in) :: exband(nbnd)
    !! Bands excluded from the calculation of overlap and projection matrices
    !
    INTEGER, INTENT (in) :: nbnd
    !! number of bands
    INTEGER, INTENT (in) :: nbndsub
    !! number of bands in the optimal subspace
    INTEGER, INTENT (in) :: nks
    !! number of kpoints in this pool
    INTEGER, INTENT (in) :: nkstot
    !! total number of kpoints
    INTEGER, INTENT (in) ::  nrr
    !! number of WS points
    INTEGER, INTENT (in) :: irvec(3, nrr)
    !! coordinates of WS points
    !
    REAL(kind=DP), INTENT (in) :: et(nbnd, nks)
    !! hamiltonian eigenvalues, coarse mesh
    REAL(kind=DP), INTENT (in) :: xk(3, nks)
    !! kpoint coordinates (cartesian in units of 2piba)
    REAL(kind=DP), INTENT (in) :: wslen(nrr)
    !! WS vectors length (alat units)
    !
    COMPLEX(kind=DP), INTENT (in) :: cu(nbnd, nbndsub, nks)
    !! rotation matrix from wannier code
    COMPLEX(kind=DP), INTENT (out) :: chw(nbndsub, nbndsub, nrr)
    !! Hamiltonian in Wannier basis
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
    INTEGER :: ibnd_m
    !! Counter on band index
    INTEGER :: i
    !! Counter on band index
    INTEGER :: nexband_tmp
    !! Number of excluded bands
    !
    REAL(kind=DP) :: rdotk
    !! $$\mathbf{r}\cdot\mathbf{k}
    REAL(kind=DP) :: tmp
    !! Maximum value of the real space Hamiltonian
    REAL(kind=dp) :: et_opt(nbnd,nks)
    !! hamiltonian eigenvalues within the outer window in the first ndimwin(ik) entries
    REAL(kind=dp) :: et_tmp(nbnd,nks)
    !! temporary array for hamiltonian eigenvalues 
    !
    COMPLEX(kind=DP) :: chs(nbndsub, nbndsub, nks)
    !! Hamiltonian in Bloch basis, coarse k-mesh
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
    ! slim down et to contain states within the outer window
    !
    nexband_tmp = 0
    DO i = 1, nbnd
      IF (exband(i)) THEN
        nexband_tmp = nexband_tmp + 1
      ENDIF  
    ENDDO
    !
    et_tmp = zero
    et_opt = zero
    IF (nexband_tmp .gt. 0) THEN 
      DO ik = 1, nks
        ibnd = 0
        DO i = 1, nbnd
          IF (exband(i)) CYCLE
          IF (lwin(i,ik)) THEN
            ibnd = ibnd + 1
            et_tmp(ibnd,ik) = et(i,ik)
          ENDIF
        ENDDO
      ENDDO
      DO ik = 1, nks
        ibnd = 0
        DO i = 1, nbnd
          IF (exband(i)) CYCLE
          ibnd = ibnd + 1
          et_opt(i,ik) = et_tmp(ibnd,ik)
        ENDDO
      ENDDO
    ELSE
      DO ik = 1, nks
        ibnd = 0
        DO i = 1, nbnd
          IF (lwin(i,ik)) THEN
            ibnd = ibnd + 1
            et_opt(ibnd,ik) = et(i,ik)
          ENDIF
        ENDDO
      ENDDO
    ENDIF
    !
    ! [Eqn. 26 of PRB 76, 165108 (2007)]
    ! H~(k) = U(k)^\dagger * H(k) * U(k)
    ! H~(k) is chs( nbndsub, nbndsub, ik )
    ! Note H~(k) is H^(W)(k) in PRB 74, 195118 (2006) notations
    !
    chs(:,:,:) = czero
    DO ik = 1, nks
      !
      DO jbnd = 1, nbndsub
        DO ibnd = 1, jbnd
          !
          ctmp = czero
          DO mbnd = 1, nbnd
            ctmp = ctmp + conjg(cu(mbnd,ibnd,ik)) * et_opt(mbnd,ik) * cu(mbnd,jbnd,ik)
          ENDDO
          !
          chs(ibnd, jbnd, ik) = ctmp
          chs(jbnd, ibnd, ik) = conjg(ctmp)
          !
        ENDDO
      ENDDO
      !
    ENDDO
    !
    CALL stop_clock ( 'Ham: step 1' )
    !
    !----------------------------------------------------------
    !  STEP 2: Fourier transform to go into Wannier basis
    !----------------------------------------------------------
    !
    ! [Eqn. 26 of PRB 76, 165108 (2007)]
    ! H(R) = (1/nk) sum_k e^{-ikR} H~(k)
    ! H(R) is chw(nbndsub, nbndsub, nrr)
    ! Note H(R) is H^(W)(R) in PRB 74, 195118 (2006) notations
    !
    CALL start_clock ( 'Ham: step 2' )
    !
    ! bring xk in crystal coordinates
    !
    CALL cryst_to_cart (nks, xk, at, -1)
    !
    chw( :, :, :) = czero
    !
    DO ir = 1, nrr
      DO ik = 1, nks
         !
         rdotk = twopi * dot_product( xk(:,ik), dble(irvec(:,ir) ))
         cfac = exp( -ci*rdotk ) / dble(nkstot)
         chw(:,:,ir) = chw(:,:,ir ) + cfac * chs(:,:,ik)
         !
      ENDDO
    ENDDO
    CALL mp_sum(chw,inter_pool_comm)
    !
    ! bring xk back into cart coord
    !
    CALL cryst_to_cart (nks, xk, bg, 1)
    !
    ! check spatial decay of Hamiltonian in Wannier basis
    !
    ! Hamiltonian matrix is in Ryd units and spatial dimensions in bohr
    ! [mind when comparing with wannier code (eV and angstrom units) with
    ! write_hr=.true.]
    !
    IF (mpime .eq. ionode_id) THEN
      !
      OPEN(unit=iudecayH,file='decay.H')
      WRITE(iudecayH, '(/3x,a/)') '#Spatial decay of Hamiltonian in Wannier basis'
      DO ir = 1, nrr
        !
        tmp =  maxval( abs( chw (:,:,ir)) )
        WRITE(iudecayH,*) wslen(ir) * celldm (1) * bohr2ang, tmp
        !
      ENDDO
      !
      ! RMDB
      DO ir = 1, nrr
        DO jbnd = 1, nbndsub
          DO ibnd = 1, nbndsub
             WRITE(iudecayH,'(5I5,6F12.6)') irvec(:,ir), ibnd, jbnd, &
                 chw(ibnd,jbnd,ir) * ryd2ev
          ENDDO
        ENDDO
      ENDDO
      !
      CLOSE(iudecayH)
      !
    ENDIF
    CALL mp_barrier(inter_pool_comm)
    !
    CALL stop_clock ( 'Ham: step 2' )
    !
    END SUBROUTINE hambloch2wan 
    !
    !--------------------------------------------------------------------------
    SUBROUTINE dmebloch2wan ( nbnd, nbndsub, nks, nkstot, dmec, xk, cu, &
       nrr, irvec, wslen, lwin, exband )
    !--------------------------------------------------------------------------
    !!
    !!  From the Dipole in Bloch representationi (coarse mesh), 
    !!  find the corresponding Dipole in Wannier representation 
    !!
    !
    USE kinds,         ONLY : DP
    USE cell_base,     ONLY : at, bg, celldm
    USE elph2,         ONLY : cdmew
    USE io_epw,        ONLY : iudecayP
    USE constants_epw, ONLY : bohr2ang, twopi, ci, czero, cone
    USE io_global,     ONLY : ionode_id
    USE mp_global,     ONLY : inter_pool_comm
    USE mp_world,      ONLY : mpime
    USE mp,            ONLY : mp_barrier, mp_sum
    ! 
    implicit none
    !
    !  input variables
    !
    INTEGER, INTENT (in) :: nbnd
    !! number of bands
    INTEGER, INTENT (in) :: nbndsub
    !! number of bands in the optimal subspace
    INTEGER, INTENT (in) :: nks
    !! number of kpoints in this pool
    INTEGER, INTENT (in) :: nkstot
    !! total number of kpoints
    INTEGER, INTENT (in) :: nrr
    !! number of WS points
    INTEGER, INTENT (in) :: irvec(3, nrr)
    !! Coordinate of Wannier space points
    !
    REAL(kind=DP), INTENT (in) :: xk(3, nks)
    !! kpoint coordinates (cartesian in units of 2piba)
    REAL(kind=DP), INTENT (in) :: wslen (nrr)
    !! WS vectors length (alat units)
    !
    COMPLEX(kind=DP), INTENT (in) :: dmec(3, nbnd, nbnd, nks)
    !! Dipole matrix elements on coarse mesh
    COMPLEX(kind=DP), INTENT (in) :: cu(nbnd, nbndsub, nks)
    !! rotation matrix from wannier code
    !
    LOGICAL, INTENT(in) :: lwin(nbnd,nks)
    !! Bands at k within outer energy window
    LOGICAL, INTENT(in) :: exband(nbnd)
    !! Bands excluded from the calculation of overlap and projection matrices
    !
    ! Local variables
    !
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
    !! Counter on band index
    INTEGER :: jbnd
    !! Counter on band index
    INTEGER :: ibnd_i
    !! Counter on band index
    INTEGER :: ibnd_j
    !! Counter on band index
    INTEGER :: nexband_tmp
    !! Number of excluded bands
    !
    REAL(kind=DP) :: rdotk
    !! $$ mathbf{r}\cdot\mathbf{k} $$
    REAL(kind=DP) :: tmp
    !! Temporary variables
    !
    COMPLEX(kind=DP) :: cps(3, nbndsub, nbndsub, nks)
    !! Dipole in smooth Bloch basis, coarse mesh
    COMPLEX(kind=DP) :: cfac
    !! $$ e^{-i\mathbf{r}\cdot\mathbf{k}} $$
    COMPLEX(kind=DP) :: dmec_utmp(nbnd,nbndsub)
    !! dmec after multiplication with the Wannier rotation matrix cu.
    COMPLEX(kind=DP) :: dmec_opt(3, nbnd, nbnd, nks)
    !! dmec computed in pmn, rescaled down if skipping lwin.
    COMPLEX(kind=DP) :: dmec_tmp(3, nbnd, nbnd, nks)
    !! temporary dmec matrix
    !
    CALL start_clock ( 'Dipole: step 1' )
    !
    !--------------------------------------------------------------
    !    STEP 0: Rescale the optical matrix on the coarse grid down
    !            This is if you skip band during the Wannierization
    !--------------------------------------------------------------
    !
    nexband_tmp = 0
    DO i = 1, nbnd
      IF (exband(i)) THEN
        nexband_tmp = nexband_tmp + 1
      ENDIF  
    ENDDO
    !
    dmec_opt = czero
    dmec_tmp = czero
    IF (nexband_tmp .gt. 0) THEN
      DO ik = 1,nks
        jbnd = 0
        DO j = 1, nbnd
          IF (exband(j)) CYCLE
          IF (lwin(j,ik)) THEN
            jbnd = jbnd + 1
            ibnd = 0
            DO i = 1, nbnd
              IF (exband(i)) CYCLE
              IF (lwin(i,ik)) THEN
                ibnd = ibnd + 1
                dmec_tmp(:,ibnd,jbnd,ik) = dmec(:,i,j,ik)
              ENDIF
            ENDDO
          ENDIF
        ENDDO
      ENDDO
      DO ik = 1,nks
        jbnd = 0
        DO j = 1, nbnd
          IF (exband(j)) CYCLE
          jbnd = jbnd + 1
          ibnd = 0
          DO i = 1, nbnd
            IF (exband(i)) CYCLE
            ibnd = ibnd + 1
            dmec_opt(:,i,j,ik) = dmec_tmp(:,ibnd,jbnd,ik)
          ENDDO
        ENDDO
      ENDDO
    ELSE
      DO ik = 1,nks
        jbnd = 0
        DO j = 1, nbnd
          IF (lwin(j,ik)) THEN
            jbnd = jbnd + 1
            ibnd = 0
            DO i = 1, nbnd
              IF (lwin(i,ik)) THEN
                ibnd = ibnd + 1
                dmec_opt(:,ibnd,jbnd,ik) = dmec(:,i,j,ik)
              ENDIF
            ENDDO
          ENDIF
        ENDDO
      ENDDO
    ENDIF
    ! 
    !----------------------------------------------------------
    !    STEP 1: rotation to optimally smooth Bloch states
    !----------------------------------------------------------
    !
    ! p~(k) = U(k)^\dagger * p(k) * U(k)
    ! p~(k) is cps( ipol, nbndsub, nbndsub, ik )
    ! Note p~(k) is p^(W)(k) in PRB 74, 195118 (2006) notations
    !
    dmec_utmp(:,:) = czero
    cps(:,:,:,:) = czero
    ! 
    DO ik = 1, nks
      DO ipol = 1, 3
        !
        ! dmec_utmp(:,:) = matmul( dmec_opt(ipol,:,:,ik), cu(:,:,ik) )
        ! cps(ipol,:,:,ik) = matmul( conjg(transpose( cu(:,:,ik))), dmec_utmp(:,:) )
        !
        CALL zgemm ('n', 'n', nbnd, nbndsub, nbnd, cone, dmec_opt(ipol,:,:,ik), &
                   nbnd, cu(:,:,ik), nbnd, czero, dmec_utmp(:,:), nbnd)
        CALL zgemm ('c', 'n', nbndsub, nbndsub, nbnd, cone, cu(:,:,ik), &
                   nbnd, dmec_utmp(:,:), nbnd, czero, cps(ipol,:,:,ik), nbndsub)
        !
      ENDDO
    ENDDO
    !
    !----------------------------------------------------------
    !  STEP 2: Fourier transform to go into Wannier basis
    !----------------------------------------------------------
    !
    ! p(R) = (1/nk) sum_k e^{-ikR} p~(k)
    ! p(R) is cdmew(ipol, nbndsub, nbndsub, nrr)
    ! Note p(R) is p^(W)(R) in PRB 74, 195118 (2006) notations
    !
    CALL start_clock ( 'Dipole: step 2' )
    !
    ! bring xk in crystal coordinates
    !
    CALL cryst_to_cart (nks, xk, at, -1)
    !
    cdmew( :, :, :, :) = czero 
    ! 
    DO ir = 1, nrr
      DO ik = 1, nks
         !
         rdotk = twopi * dot_product( xk(:,ik), dble(irvec(:,ir) ))
         cfac = exp( -ci*rdotk ) / dble(nkstot)
         cdmew(:,:,:,ir) = cdmew(:,:,:,ir) + cfac * cps(:,:,:,ik)
         !
      ENDDO
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
      OPEN(unit=iudecayP,file='decay.P')
      WRITE(iudecayP, '(/3x,a/)') '#Spatial decay of dipole in Wannier basis'
      DO ir = 1, nrr
        !
        tmp =  maxval( abs( cdmew(:, :,:,ir)) )
        WRITE(iudecayP,*) wslen(ir) * celldm(1) * bohr2ang, tmp
        !
      ENDDO
      !
      ! RMDB
      !DO ir = 1, nrr
      !  DO jbnd = 1, nbndsub
      !    DO ibnd = 1, nbndsub
      !     WRITE(iudecayP,'(5I5,6F12.6)') irvec(:,ir), ibnd, jbnd, cdmew(:,ibnd,jbnd,ir) 
      !    ENDDO
      !  ENDDO
      !ENDDO
      !     
      CLOSE(iudecayP)
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
    USE io_epw,        ONLY : iudecaydyn
    USE constants_epw, ONLY : bohr2ang, twopi, ci, czero
    USE io_global,     ONLY : ionode_id
    USE mp_world,      ONLY : mpime
    USE mp,            ONLY : mp_barrier
    USE mp_global,     ONLY : inter_pool_comm
    ! 
    implicit none
    !
    !  input variables
    !
    INTEGER, INTENT (in) :: nmodes
    !! number of branches
    INTEGER, INTENT (in) :: nq
    !! number of qpoints
    INTEGER, INTENT (in) :: nrr
    !! number of WS points
    INTEGER, INTENT (in) :: irvec(3, nrr)
    !! coordinates of WS points
    !
    REAL(kind=DP), INTENT (in) :: xk(3, nq)
    !! kpoint coordinates (cartesian in units of 2piba)
    REAL(kind=DP), INTENT (in) :: wslen(nrr)
    !! WS vectors length (alat units)
    !
    COMPLEX(kind=DP), INTENT (in) :: dynq(nmodes, nmodes, nq)
    !! dynamical matrix in bloch representation (Cartesian coordinates)
    !
    ! work variables
    !
    INTEGER :: ik
    !! Counter on k-point
    INTEGER :: ir
    !! Counter on WS points
    !
    REAL(kind=DP) :: rdotk
    !! $$ mathbf{r}\cdot\mathbf{k} $$
    REAL(kind=DP) :: tmp
    !! Temporary variables
    !
    COMPLEX(kind=DP) :: cfac
    !! $$ e^{-i\mathbf{r}\cdot\mathbf{k}} $$
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
    ! [Eqn. 29 of PRB 76, 165108 (2007)]
    ! D(R) = (1/nk) sum_k e^{-ikR} D~(k)
    ! D(R) is rdw(nmodes, nmodes, ir)
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
         rdotk = twopi * dot_product( xk( :, ik), dble(irvec( :, ir) ))
         cfac = exp( -ci*rdotk ) / dble(nq)
      !DBSP - real was commented
         rdw( :, :, ir ) = rdw( :, :, ir ) +  cfac * dynq( :, :, ik ) 
        ! rdw( :, :, ir ) = rdw( :, :, ir ) + real ( cfac * dynq( :, :, ik ) ) 
         !                                    ^^^^
         !                                 note this
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
      OPEN(unit=iudecaydyn,file='decay.dynmat')
      WRITE(iudecaydyn, '(/3x,a/)') '#Spatial decay of Dynamical matrix in Wannier basis'
      DO ir = 1, nrr
        !
        tmp =  maxval ( abs( rdw(:,:,ir)) )
        WRITE(iudecaydyn, *) wslen(ir) * celldm(1) * bohr2ang, tmp
        !
      ENDDO
      CLOSE(iudecaydyn)
    ENDIF
    CALL mp_barrier(inter_pool_comm)
    !
    END SUBROUTINE dynbloch2wan
    !-----------------------------------------------------
    !--------------------------------------------------------------------------
    SUBROUTINE vmebloch2wan ( nbnd, nbndsub, nks, nkstot, xk, cu, &
       nrr, irvec, wslen, lwin, exband )
    !--------------------------------------------------------------------------
    !!
    !!  Calculate the velocity matrix elements in the Wannier basis
    !!  at no point do we actually have the coarse mesh v-ME. 
    !!
    !! RM 03/2018: debugged and updated
    !
    !--------------------------------------------------------------------------
    !
    USE kinds,     ONLY : DP
    USE cell_base, ONLY : at, bg, celldm
    USE elph2,     ONLY : cvmew
    USE constants_epw, ONLY : twopi, one, zero, ci, czero, cone, bohr2ang
    USE io_epw,    ONLY : iummn, iubvec, iudecayv
    USE io_files,  ONLY : prefix
    USE io_global, ONLY : ionode_id, stdout
    USE mp_global, ONLY : inter_pool_comm, my_pool_id
    USE mp,        ONLY : mp_barrier, mp_sum
    USE mp_world,  ONLY : mpime
    implicit none
    !
    !  input variables
    !
    INTEGER, INTENT (in) :: nbnd
    !! number of bands
    INTEGER, INTENT (in) :: nbndsub
    !! number of bands in the optimal subspace
    INTEGER, INTENT (in) :: nks
    !! number of kpoints in this pool
    INTEGER, INTENT (in) :: nkstot
    !! total number of kpoints
    INTEGER, INTENT (in) :: nrr
    !! number of WS points
    INTEGER, INTENT (in) :: irvec(3, nrr)
    !! Coordinate of Wannier space points
    !
    REAL(kind=DP), INTENT (in) :: xk(3, nks)
    !! kpoint coordinates (cartesian in units of 2piba)
    REAL(kind=DP), INTENT (in) :: wslen(nrr)
    !! WS vectors length (alat units)
    !
    COMPLEX(kind=DP), INTENT (in) :: cu(nbnd, nbndsub, nks)
    !! rotation matrix from wannier code
    !
    LOGICAL, INTENT(in) :: lwin(nbnd,nks)
    !! Bands at k within outer energy window
    LOGICAL, INTENT(in) :: exband(nbnd)
    !! Bands excluded from the calculation of overlap and projection matrices
    !
    ! Local variables
    !
    INTEGER :: ipol
    !! Counter on polarization
    INTEGER :: ik
    !! Counter on k-point
    INTEGER :: ir
    !! Counter on WS points
    INTEGER :: ikstart
    !! Starting ik for this pool
    INTEGER :: ikstop
    !! Ending ik for this pool
    INTEGER :: ibnd
    !! Counter on band index
    INTEGER :: jbnd
    !! Counter on band index
    INTEGER :: i
    !! Counter on band index
    INTEGER :: j
    !! Counter on band index
    INTEGER :: ibnd_i
    !! Counter on band index
    INTEGER :: ibnd_j
    !! Counter on band index
    INTEGER :: nnb
    !! total number of neighbours for each k-point
    INTEGER :: ib
    !! Counter on b-vectors
    INTEGER :: ipool
    !! Current pool
    INTEGER :: nkb
    !! in-pool index of k+b-point
    INTEGER :: nkb_abs
    !! absolute index of k+b-point
    INTEGER :: nkk
    !! in-pool index of k-point
    INTEGER :: nkk_abs
    !! absolute index of k-point
    INTEGER :: ios
    !! Integer variable for I/O control
    INTEGER :: nexband_tmp
    !! Number of excluded bands
    !
    REAL(kind=DP), ALLOCATABLE :: bvec(:,:,:)
    !! b-vectors connecting each k-point to its nearest neighbors
    REAL(kind=DP), ALLOCATABLE :: wb(:)
    !! weight of the nnb-th b-vector (ordered along shells)
    REAL(kind=DP) :: rdotk
    !! $$ mathbf{r}\cdot\mathbf{k} $$
    REAL(kind=DP) :: tmp
    !! Temporary variables
    COMPLEX(kind=DP) :: cfac
    !! $$ e^{-i\mathbf{r}\cdot\mathbf{k}} $$
    COMPLEX(kind=DP) :: cu_big (nbnd, nbndsub, nkstot)
    !! rotation matrix from wannier code
    REAL(kind=DP) ::  b_tmp(3)
    !! temporary b-vectors
    REAL(kind=DP) :: zero_vect(3)
    !! temporary zero vector
    REAL(kind=DP) :: delta
    !! \delta_nm = 1 if n .eq. m and 0 if n .neq. m
    !
    COMPLEX(kind=DP) :: Apos(3,nbndsub,nbndsub,nks)
    !! A^W_{mn,\alpha}(k)
    COMPLEX(kind=DP), ALLOCATABLE :: M_mn(:,:,:,:)
    !! M_mn(k,b)
    COMPLEX(kind=DP), ALLOCATABLE :: m_mat_opt(:,:,:,:)
    !! M_mn(k,b) computed in mmn, rescaled down if skipping lwin.
    COMPLEX(kind=DP), ALLOCATABLE :: m_mat_tmp(:,:,:,:)
    !! temporary M_mn matrix
    COMPLEX(kind=DP), ALLOCATABLE :: cvs(:,:,:,:)
    !! M_mn in smooth Bloch basis, coarse k-mesh
    COMPLEX(kind=DP) :: M_mn_utmp(nbnd,nbndsub)
    !! M_mn after multiplication with the Wannier rotation matrix cu.
    COMPLEX(kind=DP) :: ctmp
    !! Temporary variable to store M_mn
    !
    CHARACTER (len=256) :: tempfile
    !! Temporary file
    !
    ! setup rotation matrix - we need access to all for the k+b
    cu_big = czero
    CALL ckbounds(ikstart, ikstop)
    cu_big(:,:,ikstart:ikstop) = cu(:,:,:)
    CALL mp_sum(cu_big,inter_pool_comm)
    !
    !--------------------------------------------------------------
    !  STEP 0: Read in b-vectors bvec and their weights wb, and
    !          M_mn matrix
    !--------------------------------------------------------------
    !
    ! RM - bvec can be writen on file by making a small change in
    ! W90/hamiltonian.F90/hamilotonian_write_rmn
    !
    tempfile=trim(prefix)//'.bvec'
    OPEN(iubvec, file=tempfile, action='read', iostat=ios)
    IF (ios /= 0) THEN
      !
      ! if it doesn't exist, then we just set the bvec and wb to zero
      !
      nnb = 1
      ALLOCATE ( bvec(3,nnb,nkstot), wb(nnb) )
      bvec = zero
      wb   = zero
    ELSE
      READ(iubvec,*) nnb
      ALLOCATE ( bvec(3,nnb,nkstot), wb(nnb) )
      DO ik = 1, nkstot
        DO ib = 1, nnb
          READ(iubvec,*) bvec(:,ib,ik), wb(ib)
        ENDDO
      ENDDO
      CLOSE(iubvec)
    ENDIF
    !
    ! b-vectors in units of Ang^-1 and wb-weights in units of Ang^2
    ! when read from file
    bvec = bvec * bohr2ang
    wb = wb / bohr2ang**2
    !
    !  read Mmn for velocity calculation
    !
    ! RM - M_mn matrix is writen on file in pw2wan90epw.f90/compute_mmn_para
    !    - dimensions of M_mn are M_mn(nbnd, nbnd, nnb, nkstot)
    !
    ALLOCATE ( M_mn(nbnd, nbnd, nnb, nkstot) )
    M_mn = czero
    !
    IF (mpime.eq.ionode_id) THEN
      tempfile=trim(prefix)//'.mmn'
      OPEN(iummn, file=tempfile, status = 'old', form = 'formatted', iostat=ios)
      !
      IF (ios /= 0) THEN
        ! if it doesn't exist, then we just set the mmn to zero
        CALL errore ('vmebloch2wan','error opening' // tempfile, 0)
      ELSE
        !
        DO ik = 1, nkstot
          DO ib = 1, nnb
            DO jbnd = 1, nbnd
              DO ibnd = 1, nbnd
                READ(iummn,*) M_mn(ibnd,jbnd,ib,ik)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
        CLOSE(iummn)
        !
      ENDIF
    ENDIF
    CALL mp_sum(M_mn,inter_pool_comm)
    !
    CALL start_clock ( 'Velocity: step 1' )
    !
    ! slim down M_mn matrix to contain states within the outer window
    !
    nexband_tmp = 0
    DO i = 1, nbnd
      IF (exband(i)) THEN
        nexband_tmp = nexband_tmp + 1
      ENDIF  
    ENDDO
    !
    ALLOCATE ( m_mat_opt(nbnd, nbnd, nnb, nks) )
    ALLOCATE ( m_mat_tmp(nbnd, nbnd, nnb, nks) )
    m_mat_opt(:,:,:,:) = czero
    m_mat_tmp(:,:,:,:) = czero
    zero_vect(:) = zero
    !
    IF (nexband_tmp .gt. 0) THEN
      DO ik = 1, nks
        CALL ktokpmq ( xk(:,ik), zero_vect, +1, ipool, nkk, nkk_abs)
        !
        jbnd = 0
        DO j = 1, nbnd
          IF (exband(j)) CYCLE
          IF (lwin(j,ik)) THEN
            jbnd = jbnd + 1
            ibnd = 0
            DO i = 1, nbnd
              IF (exband(i)) CYCLE
              IF (lwin(i,ik)) THEN
                ibnd = ibnd + 1
                m_mat_tmp(ibnd,jbnd,:,ik) = M_mn(i,j,:,nkk_abs)
              ENDIF
            ENDDO
          ENDIF
        ENDDO
      ENDDO
      DO ik = 1,nks
        jbnd = 0
        DO j = 1, nbnd
          IF (exband(j)) CYCLE
          jbnd = jbnd + 1
          ibnd = 0
          DO i = 1, nbnd
            IF (exband(i)) CYCLE
            ibnd = ibnd + 1
            m_mat_opt(i,j,:,ik) = m_mat_tmp(ibnd,jbnd,:,ik)
          ENDDO
        ENDDO
      ENDDO
    ELSE
      DO ik = 1, nks
        CALL ktokpmq ( xk(:,ik), zero_vect, +1, ipool, nkk, nkk_abs)
        !
        jbnd = 0
        DO j = 1, nbnd
          IF (lwin(j,ik)) THEN
            jbnd = jbnd + 1
            ibnd = 0
            DO i = 1, nbnd
              IF (lwin(i,ik)) THEN
                ibnd = ibnd + 1
                m_mat_opt(ibnd,jbnd,:,ik) = M_mn(i,j,:,nkk_abs)
              ENDIF
            ENDDO
          ENDIF
        ENDDO
      ENDDO
    ENDIF
    !
    DEALLOCATE (M_mn)
    DEALLOCATE (m_mat_tmp)
    !
    !----------------------------------------------------------
    ! STEP 1: Calculate A^(W)_{mn,\alpha}(k) [Eqn. 44 of PRB 74, 195118 (2006)]
    ! A^(W)_{mn,\alpha}(k) = i \sum_b wb b_{\alpha} (<u^(W)_nk | u^(W)_m(k+b)> - delta_mn)
    !----------------------------------------------------------
    !
    ! <u^(W)_nk | u^(W)_m(k+b)> = U(k)^\dagger * M_mn(k) * U(k+b)
    !
    ! <u^(W)_nk | u^(W)_m(k+b)> is cvs(nbndsub,nbndsub,nnb,nks)
    !
    ALLOCATE ( cvs(nbndsub,nbndsub,nnb,nks) )
    !
    cvs(:,:,:,:) = czero
    M_mn_utmp(:,:) = czero
    !
    DO ik = 1, nks
      !
      CALL ktokpmq ( xk(:,ik), zero_vect, +1, ipool, nkk, nkk_abs)
      !
      DO ib = 1, nnb
        !
        ! bring bvec to units of 2piba since xk is cartesian units of 2piba
        b_tmp(:) = celldm(1) / (twopi) * bvec(:,ib,nkk_abs)
        CALL ktokpmq ( xk(:,ik), b_tmp(:), +1, ipool, nkb, nkb_abs)
        !
        ! M_mn_utmp(:,:) = matmul( m_mat_opt(:,:,ib,ik), cu_big(:,:,nkb_abs) )
        ! cvs(:,:,ib,ik) = matmul( conjg(transpose(cu(:,:,ik))), M_mn_utmp(:,:) )
        !
        CALL zgemm ('n', 'n', nbnd, nbndsub, nbnd, cone, m_mat_opt(:,:,ib,ik), &
                   nbnd, cu_big(:,:,nkb_abs), nbnd, czero, M_mn_utmp(:,:), nbnd)
        CALL zgemm ('c', 'n', nbndsub, nbndsub, nbnd, cone, cu(:,:,ik), &
                   nbnd, M_mn_utmp(:,:), nbnd, czero, cvs(:,:,ib,ik), nbndsub)
      ENDDO
      !
    ENDDO
    !
    DEALLOCATE ( m_mat_opt )
    !
    ! A^(W)_{mn,\alpha}(k) is Apos(3,nbndsub,nbndsub,nks)
    !
    Apos = czero
    !
    DO ik = 1, nks
      !
      CALL ktokpmq ( xk(:,ik), zero_vect, +1, ipool, nkk, nkk_abs)
      !
      DO ib = 1, nnb
        !
        DO jbnd = 1, nbndsub
          DO ibnd = 1, nbndsub
            !
            delta = zero
            IF ( jbnd == ibnd ) delta = one
            DO ipol = 1, 3
              Apos(ipol,ibnd,jbnd,ik) = Apos(ipol,ibnd,jbnd,ik) + &
                 ci * wb(ib) * bvec(ipol,ib,nkk_abs) * ( cvs(ibnd,jbnd,ib,ik) - delta )
            ENDDO
            !
          ENDDO
        ENDDO
        !
      ENDDO
      !
    ENDDO
    !
    DEALLOCATE ( bvec )
    DEALLOCATE ( wb )
    !
    CALL stop_clock ( 'Velocity: step 1' )
    !
    !----------------------------------------------------------
    ! STEP 2: Fourier transform to go into Wannier basis
    !----------------------------------------------------------
    !
    ! Calculate position matrix [Eqn. 43 of PRB 74, 195118 (2006)]
    ! r_{\alpha}(R) = (1/nk) \sum_k e^(-ikR) A^(W)_{mn,\alpha}(k)
    !
    CALL start_clock ( 'Velocity: step 2' )
    !
    ! bring xk in crystal coordinates
    !
    CALL cryst_to_cart (nks, xk, at, -1)
    !
    ! r_{\alpha}(R) is cvmew(3,nbndsub,nbndsub,nrr)
    !
    cvmew(:,:,:,:) = czero
    ! 
    DO ir = 1, nrr
      DO ik = 1, nks
         !
         rdotk = twopi * dot_product( xk(:,ik), dble(irvec(:,ir) ))
         cfac = exp( -ci*rdotk ) / dble(nkstot)
         cvmew(:,:,:,ir) = cvmew(:,:,:,ir) + cfac * Apos(:,:,:,ik)
         !
      ENDDO
    ENDDO
    !
    CALL mp_sum(cvmew,inter_pool_comm) 
    !
    ! bring xk back into cart coord
    !
    CALL cryst_to_cart (nks, xk, bg, 1)
    !
    ! check spatial decay of position matrix elements in Wannier basis
    !
    ! position matrix cvmew and spatial dimensions are in units of bohr
    ! [mind when comparing with wannier code (angstrom units) with write_rmn=.true.]
    !
    IF (mpime.eq.ionode_id) then
      OPEN(unit=iudecayv,file='decay.v')
      WRITE(iudecayv, '(/3x,a/)') '#Spatial decay of Velocity matrix element in Wannier basis'
      DO ir = 1, nrr
        !
        tmp =  maxval ( abs( cvmew(:,:,:,ir)) )
        WRITE(iudecayv, *) wslen(ir) * celldm(1) * bohr2ang, tmp
        !
      ENDDO
      !
      ! RMDB
      !DO ir = 1, nrr
      !  DO jbnd = 1, nbndsub
      !    DO ibnd = 1, nbndsub
      !      WRITE(iudecayv,'(5I5,6F12.6)') irvec(:,ir), ibnd, jbnd, cvmew(:,ibnd,jbnd,ir) * bohr2ang
      !    ENDDO
      !  ENDDO
      !ENDDO
      !
      CLOSE(iudecayv)
    ENDIF
    CALL mp_barrier(inter_pool_comm)
    !
    CALL stop_clock ( 'Velocity: step 2' )
    !
    WRITE(stdout,'(/5x,a)') 'Velocity matrix elements calculated'
    WRITE(stdout,'(a)') ' '
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
    !! Number of kpoints in this pool
    INTEGER, INTENT (in) :: nrr
    !! Number of WS points
    INTEGER, INTENT (in) :: nkstot
    !! Total number of kpoints
    INTEGER, INTENT (in) :: irvec(3, nrr)
    !! Coordinates of WS points
    !
    REAL(kind=DP), INTENT (in) :: xk(3, nks)
    !! kpoint coordinates (cartesian in units of 2piba)
    REAL(kind=DP), INTENT (in) :: wslen(nrr)
    !! WS vectors length (alat units)
    !
    COMPLEX(kind=DP), INTENT (in) :: cu(nbnd, nbndsub, nks)
    !! Rotation matrix from wannier code
    COMPLEX(kind=DP), INTENT (in) :: cuq(nbnd, nbndsub, nks)
    !! Rotation matrix from wannier code
    COMPLEX(kind=DP), INTENT (in) :: epmatk( nbnd, nbnd, nks)
    !! e-p matrix in bloch representation, coarse mesh
    !
    ! output variables
    !
    COMPLEX(kind=DP), INTENT(out) :: epmatw( nbndsub, nbndsub, nrr)
    !!  e-p matrix  in Wannier basis
    !
    ! Work variables
    !
    INTEGER :: ik
    !! Counter on k-point
    INTEGER :: ir
    !! Counter on WS points
    INTEGER :: ibnd, jbnd
    !
    REAL(kind=DP) :: rdotk
    !! $$ mathbf{r}\cdot\mathbf{k} $$
    REAL(kind=DP) :: tmp
    !! Temporary variables
    !
    COMPLEX(kind=DP) :: cfac
    !! $$ e^{-i\mathbf{r}\cdot\mathbf{k}} $$
    COMPLEX(kind=DP) :: epmats(nbndsub, nbndsub, nks)
    !!  e-p matrix  in smooth Bloch basis, coarse mesh
    COMPLEX(kind=DP) :: eptmp(nbndsub, nbnd)
    !!  e-p matrix, temporary
    !
    !----------------------------------------------------------
    !  STEP 1: rotation to optimally smooth Bloch states
    !----------------------------------------------------------
    !
    ! [Eqn. 24 of PRB 76, 165108 (2007)]
    ! g~(k,q) = U(k+q)^\dagger * g(k,q) * U(k)
    !
    ! g(k,q)   is epmatk (ibnd, jbnd, ik)
    ! g~(k,q)  is epmats (ibnd, jbnd, ik)
    !
    CALL start_clock ( 'ep: step 1' )
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
    ! [Eqn. 24 of PRB 76, 165108 (2007)]
    ! g(R_e,q) = (1/nkc) sum_k e^{-ikR_e} g~(k,q)
    ! g(R_e,q) is epmatw (nbndsub,nbndsub,ir)
    !
    CALL start_clock ( 'ep: step 2' )
    !
    epmatw(:, :, :) = czero
    !
    ! bring xk in crystal coordinates
    !
    CALL cryst_to_cart (nks, xk, at, -1)
    !
    DO ir = 1, nrr
       DO ik = 1, nks
         !
         rdotk = twopi * dot_product( xk( :, ik), dble(irvec( :, ir) ))
         cfac = exp( -ci*rdotk ) / dble(nkstot)
         epmatw( :, :, ir) = epmatw( :, :, ir) + cfac * epmats( :, :, ik)
         !
       ENDDO
    ENDDO
    !
    CALL mp_sum(epmatw,inter_pool_comm)  
    !
    ! bring xk back into cart coord
    !
    CALL cryst_to_cart (nks, xk, bg, 1)
    !
    !
    !  Check spatial decay of matrix elements in Wannier basis
    !  the unit in r-space is angstrom, and I am plotting 
    !  the matrix for the first mode only
    !
    IF (mpime.eq.ionode_id) THEN
      OPEN(unit=iuwane,file='decay.epwane')
      WRITE(iuwane, '(a)') '# Spatial decay of e-p matrix elements in Wannier basis'
      DO ir = 1, nrr
        ! 
        tmp =  maxval ( abs(epmatw(:,:,ir)) ) 
        WRITE(iuwane, *) wslen(ir) * celldm(1) * bohr2ang, tmp
        !
      ENDDO
      !
      ! RMDB
      !DO ir = 1, nrr
      !  DO jbnd = 1, nbndsub
      !    DO ibnd = 1, nbndsub
      !      WRITE(iuwane,'(5I5,2F12.6)') irvec(:,ir), ibnd, jbnd, epmatw(jbnd,ibnd,ir)
      !    ENDDO
      !  ENDDO
      !ENDDO
      !
      CLOSE(iuwane)
    ENDIF
    !
    CALL stop_clock ( 'ep: step 2' )
    !
    END SUBROUTINE ephbloch2wane  
    ! 
    !--------------------------------------------------------------------------
    SUBROUTINE ephbloch2wanp ( nbnd, nmodes, xk, nq, irvec_k, irvec_g, &
      nrr_k, nrr_g, epmatwe )
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
    ! 
    implicit none
    !
    !  Input variables 
    !
    INTEGER, INTENT(in) :: nbnd
    !! Number of electronic bands
    INTEGER, INTENT(in) :: nmodes
    !! number of branches
    INTEGER, INTENT(in) :: nq 
    !! number of qpoints
    INTEGER, INTENT(in) :: nrr_k
    !! number of electronic WS points
    INTEGER, INTENT(in) :: nrr_g
    !! number of el-ph WS points
    INTEGER, INTENT(in) :: irvec_k(3, nrr_k)
    !! Coordinates of real space vector for electrons
    INTEGER, INTENT(in) :: irvec_g(3, nrr_g)
    !! Coordinates of real space vector for electron-phonon
    !
    REAL(kind=DP), INTENT(in) :: xk(3, nq)
    !! Kpoint coordinates (cartesian in units of 2piba) 
    ! 
    COMPLEX(kind=DP), INTENT(in) :: epmatwe(nbnd, nbnd, nrr_k, nmodes, nq)
    !! EP matrix in electron-wannier representation and phonon bloch representation
    !!   (Cartesian coordinates)
    !
    ! EP matrix in electron-wannier representation and phonon-Wannier  representation
    !
    ! Work variables
    !
    INTEGER :: iq
    !! Counter on q-point
    INTEGER :: ir
    !! Counter on WS points
    INTEGER :: ire
    !! Counter on WS points
    !
    REAL(kind=DP) :: rdotk
    !! $$ mathbf{r}\cdot\mathbf{k} $$
    REAL(kind=DP) :: tmp
    !! Temporary variables
    !
    REAL(kind=DP) :: rvec1(3), rvec2(3), len1, len2
    !
    COMPLEX(kind=DP) :: cfac
    !! $$ e^{-i\mathbf{r}\cdot\mathbf{k}} $$
    !
    !----------------------------------------------------------
    !  Fourier transform to go into Wannier basis
    !----------------------------------------------------------
    !
    ! [Eqn. 24 of PRB 76, 165108 (2007)]
    ! g(R_e,R_p) = (1/nq) sum_q e^{-iqR_p} g(R_e,q)
    !
    ! bring xk in crystal coordinates
    !
    CALL cryst_to_cart (nq, xk, at, -1)
    !
    epmatwp = czero
    ! 
    DO ir = 1, nrr_g
      !
      DO iq = 1, nq
        !
        rdotk = twopi * dot_product( xk( :, iq), dble(irvec_g( :, ir) ))
        cfac = exp( -ci*rdotk ) / dble(nq)
        epmatwp(:,:,:,:,ir) = epmatwp(:,:,:,:,ir) + cfac * epmatwe(:,:,:,:,iq)
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
        DO ire = 1, nrr_k
          !
          rvec1 = dble(irvec_k(1,ire))*at(:,1) + &
                  dble(irvec_k(2,ire))*at(:,2) + &
                  dble(irvec_k(3,ire))*at(:,3)
          rvec2 = dble(irvec_g(1,ir))*at(:,1) + &
                  dble(irvec_g(2,ir))*at(:,2) + &
                  dble(irvec_g(3,ir))*at(:,3)
          len1 = sqrt(rvec1(1)**2.d0+rvec1(2)**2.d0+rvec1(3)**2.d0)
          len2 = sqrt(rvec2(1)**2.d0+rvec2(2)**2.d0+rvec2(3)**2.d0)
          tmp =  maxval ( abs( epmatwp (:, :, ire, :, ir) ) )
          !
          ! rvec1 : electron-electron0 distance
          ! rvec2 : phonon - electron0 distance
          !
          WRITE(iuwanep, '(5f15.10)') len1 * celldm(1) * bohr2ang, &
                                  len2 * celldm(1) * bohr2ang, tmp
        ENDDO
        IF (ir.eq.nrr_g) CLOSE(iuwanep)
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
    SUBROUTINE ephbloch2wanp_mem ( nbnd, nmodes, xk, nq, irvec_k, irvec_g, &
      nrr_k, nrr_g, epmatwe )
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
    INTEGER, INTENT(in) :: nrr_k
    !! number of electronic WS points
    INTEGER, INTENT(in) :: nrr_g
    !! number of el-h WS points
    INTEGER, INTENT(in) :: nmodes
    !! number of branches
    INTEGER, INTENT(in) :: nq
    !! number of qpoints
    INTEGER, INTENT(in) :: irvec_k (3, nrr_k)
    !! Coordinates of real space vector 
    INTEGER, INTENT(in) :: irvec_g (3, nrr_g)
    !! Coordinates of real space vector 
    !
    REAL(kind=DP), INTENT(in) :: xk(3, nq)
    !! Kpoint coordinates (cartesian in units of 2piba) 
    ! 
    COMPLEX(kind=DP), INTENT(in) :: epmatwe (nbnd, nbnd, nrr_k, nmodes)
    !! EP matrix in electron-wannier representation and phonon bloch representation
    !!   (Cartesian coordinates)
    !
    ! EP matrix in electron-wannier representation and phonon-Wannier representation
    !
    ! work variables
    !
    INTEGER :: iq
    !! Counter on q-point
    INTEGER :: ir
    !! Counter on WS points
    INTEGER :: ire
    !! Counter on WS points
    !
    REAL(kind=DP) :: rdotk
    !! $$ mathbf{r}\cdot\mathbf{k} $$
    REAL(kind=DP) :: tmp
    !! Temporary variables
    !
    REAL(kind=DP) :: rvec1(3), rvec2(3), len1, len2
    !
    COMPLEX(kind=DP) :: cfac
    !! $$ e^{-i\mathbf{r}\cdot\mathbf{k}} $$
    !
    COMPLEX(KIND=DP), ALLOCATABLE :: epmatwp_mem(:,:,:,:)
    !!  e-p matrix in Wannier basis
    !
    ALLOCATE (epmatwp_mem( nbnd, nbnd, nrr_k, nmodes))
    ! 
    !----------------------------------------------------------
    !  Fourier transform to go into Wannier basis
    !----------------------------------------------------------
    !
    ! [Eqn. 24 of PRB 76, 165108 (2007)]
    ! g(R_e,R_p) = (1/nq) sum_q e^{-iqR_p} g(R_e,q)
    !
    ! bring xk in crystal coordinates
    !
    CALL cryst_to_cart (nq, xk, at, -1)
    !
    DO ir = 1, nrr_g
      !
      epmatwp_mem = czero
      ! 
      DO iq = 1, nq
         !
         ! direct read of epmatwe for this iq 
         CALL rwepmatw ( epmatwe, nbnd, nrr_k, nmodes, iq, iunepmatwe, -1)
         !
         rdotk = twopi * dot_product( xk( :, iq), dble(irvec_g( :, ir) ))
         cfac = exp( -ci*rdotk ) / dble(nq)
         epmatwp_mem = epmatwp_mem + cfac * epmatwe
         !
      ENDDO
      !
      ! direct write of epmatwp_mem for this ir 
      CALL rwepmatw (epmatwp_mem, nbnd, nrr_k, nmodes, ir, iunepmatwp, +1)
      !  check spatial decay of e-p matrix elements in wannier basis - electrons
      !  + phonons
      !
      !  we plot: R_e, R_p, max_{m,n,nu} |g(m,n,nu;R_e,R_p)|
      !
      IF (mpime == ionode_id) THEN
        IF (ir == 1) OPEN(unit=iuwanep, file='decay.epmat_wanep', status='unknown')
        IF (ir == 1) WRITE(iuwanep, '(a)') '#  R_e,    R_p, max_{m,n,nu} |g(m,n,nu;R_e,R_p)| '
        DO ire = 1, nrr_k
          !
          rvec1 = dble(irvec_k(1,ire))*at(:,1) + &
                  dble(irvec_k(2,ire))*at(:,2) + &
                  dble(irvec_k(3,ire))*at(:,3)
          rvec2 = dble(irvec_g(1,ir))*at(:,1) + &
                  dble(irvec_g(2,ir))*at(:,2) + &
                  dble(irvec_g(3,ir))*at(:,3)
          len1 = sqrt(rvec1(1)**2.d0+rvec1(2)**2.d0+rvec1(3)**2.d0)
          len2 = sqrt(rvec2(1)**2.d0+rvec2(2)**2.d0+rvec2(3)**2.d0)
          tmp =  maxval ( abs( epmatwp_mem(:, :, ire, :) ) )
          !
          ! rvec1 : electron-electron0 distance
          ! rvec2 : phonon - electron0 distance
          !
          WRITE(iuwanep, '(5f15.10)') len1 * celldm(1) * bohr2ang, &
                                  len2 * celldm(1) * bohr2ang, tmp
        ENDDO
        IF (ir == nrr_g) CLOSE(iuwanep)
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
