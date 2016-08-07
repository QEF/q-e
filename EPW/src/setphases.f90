  !                                                                           
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino 
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  !-----------------------------------------------------------------------
  subroutine setphases ( kunit, ik0, nng, unimat )
  !---------------------------------------------------------------------
  !
  ! This subroutine gives the rotations which bring the
  ! eigenstates of the Hamiltonian (evc) into a uniquely defined
  ! gauge (independent of the machine or the numerical libraries)
  !
  ! This is done by diagonalizing a fake perturbation within the
  ! degenerate hamiltonian manifold. The eigenstates are assumed
  ! to be labeled in order of increasing energy (as it is usually 
  ! the case).
  !  
  ! The first step is the diagonalization of the perturbation.
  ! This still leaves a phase arbitrariness, which is the fixed
  ! by requiring some c(G) be real.
  !
  ! It should be :  |psi^\prime_i>  =  \sum_j U_{ji} * |psi_j>
  ! and A^\prime = U^\dagger * A * U = matmul(conjg(transpose(U)),matmul(A,U)) 
  ! with U = unimat
  !
  ! Only for 1 proc/pool (This limitation can be removed)
  !
  ! ----------------------------------------------------------------------
  !
  USE mp_global,            ONLY : my_pool_id, nproc_pool, &
                                   me_pool, intra_pool_comm
  USE mp,                   ONLY : mp_barrier, mp_sum
  USE io_global,            ONLY : stdout
  USE kinds,                ONLY : DP
  USE wvfct,                ONLY : et
  USE phcom,                ONLY : evq 
  USE qpoint,               ONLY : igkq
  USE control_flags,        ONLY : iverbosity
  USE pwcom,                ONLY : igk 
  USE constants_epw,        ONLY : ci, czero
  USE wavefunctions_module, ONLY : evc
  USE fft_base,             ONLY : dfftp, dffts
  USE wvfct,                ONLY : nbnd   
  USE gvecs,                ONLY : nls
  USE fft_interfaces,       ONLY : fwfft, invfft

!$$
  use cell_base,            ONLY : omega, pi
  ! necessary in generating the fake perturbation
!$$
  !
  implicit none
  !
  INTEGER, PARAMETER :: nnglen = 100, ndig = 5
  !   reduced size of evc, evq, this is just for setting individual phases 
  !     this may be taken = 1 in principle, but if C(G1) = 0 we are in trouble, so...
  !   number of digits used in integer comparisons to decide the largest c(G) in the set
  COMPLEX(KIND=DP), POINTER :: evtmp(:,:)     
  !   pointer to wavefunctions in the PW basis (either evc or evq)
  INTEGER, POINTER :: igktmp(:)
  !   correspondence k+G <-> G
  INTEGER :: nng, kunit, ik0, nset, ndeg(nbnd)
  !   size of evc, evq ( this is npw or npwq from elphel2) 
  !   granularity of k point distribution
  !   the current k point
  !   number of degenerate subsets for this k point
  !   degeneracy of each subset
  COMPLEX(kind=DP) :: c (nnglen,nbnd) , unimat (nbnd,nbnd), cnew (nnglen,nbnd)
  !   the chunk of evc used for setting the phases
  !   the global rotation matrix
  !   the rotated chunks
  COMPLEX(kind=DP), ALLOCATABLE ::  u (:,:)
  !   the rotation matrix for the degenarate subset
  REAL(kind=DP) :: deltav (dffts%nnr)
  !   the fake (local) perturbation in real space, it is real to guarantee Hermiticity
  REAL(kind=DP), PARAMETER :: eps = 1.d-5, epsc = 1.d-8
  !   threshold for deciding whether two states are degenerate
  !   threshold for deciding whether c(G) = 0
  !
  ! variables for lapack ZHPEVX 
  !
  integer                       :: neig, info
  integer, allocatable          :: ifail(:), iwork(:)
  real(kind=DP), allocatable    :: w(:), rwork(:)
  complex(kind=DP), allocatable :: up(:), cwork(:), cz(:,:)
  !
  ! work variables
  !
  complex(kind=DP) :: unitcheck (nbnd, nbnd), ctmp
  real(kind=DP) :: theta, tmp
  integer :: ir, ig, igmax, ibnd, jbnd, mbnd, pbnd, iset, n0, n1, n2, &
     nsmall, ibndg, jbndg, itmp, imaxnorm
  logical :: tdegen
  complex(kind=DP) :: aux1 (dffts%nnr), deltavpsi (nng)
  complex(kind=DP) :: ZDOTC

!$$ parameters used in generating deltav(ir)
  integer :: i,j,k
  integer :: na1,na2,na3,na4,na5
  ! na1 ~ na3 : initial seed for perturbation
  !
  ! na4: increment of the real space indices corresponding to
  !  the atomic scale variation, here set to 0.5 bohr. This is
  !  a very important criterion for the fake perturbation.
  !  If the variation of the fake perturbation is too rapid or
  !  too slow, it will not capture the variation of the electronic
  !  wavefunctions, which again is in atomic scale.
  !
  ! na5: integers between (0 ~ na5-1) is assigned first and
  !  then renormalized such that the eigenvalues of the fake
  !  perturbation is a number between 0 and 1.
!$$

  ! In the case of multiple procs per pool, the fake perturbation defined 
  ! by deltav (ir) = ir will depend on the division of nr1x*nr2x*nr3x
  ! (This can be removed by defined a serious perturbation below...)
  !
  IF (nproc_pool>1) call errore &
       ('setphases', 'only one proc per pool to guarantee the same gauge', 1)
  !
  ! initialize
  !
  ! if we are working with a k (k+q) point, evtmp points to evc (evq)
  ! 
  IF ( (kunit.eq.1) .or. ((kunit.eq.2).and.(ik0-2*(ik0/2).eq.1)) ) then
    evtmp  => evc
    igktmp => igk
  ELSE
    evtmp  => evq
    igktmp => igkq
  ENDIF
  c = evtmp(1:nnglen,:)
  !  
  ! build the fake perturbation
  !
  ! This is crap but it works fine. I think it is important to normalize 
  ! the pert to 1 so as to easily keep control on the eigenvalues.
  ! According to wikipedia the word crap is from latin "crappa".
  ! If you do something better you should take into account the 
  ! parallalization on nrxxs.
  !  
!  do ir = 1, nrxxs
!    deltav(ir) = dble(ir)/dble(nrxxs)
!  enddo
  !
  ! the above is not ok for my BC53 structure... A better choice:
  ! read dvscf (without the bare term) for a few patterns.
  ! The spin index is irrelevant and is kept only for compatibility with davcio_drho.
  ! To get a unique ordering independent of q, the same deltav must be used!
  ! Therefore in input we supply fildvscf0 (the fildvscf calculated, say, at gamma)
  !
!$$  call davcio_drho ( v1,  lrdrho, iudvscf0,       1, -1 )
!$$  call davcio_drho ( v2,  lrdrho, iudvscf0, 3*nat/2, -1 )
!$$  call davcio_drho ( v3,  lrdrho, iudvscf0, 3*nat  , -1 )
!$$  deltav = real ( v1(:,1) + v2(:,1) + v3(:,1) )
!$$  deltav = deltav ** 3.d0

  !
  nset = 1 
  DO ibnd = 1, nbnd
    ndeg (ibnd) = 1
    DO jbnd = 1, nbnd
       unimat (ibnd, jbnd) = czero
    ENDDO
    unimat (ibnd, ibnd) = 1.d0 
  ENDDO
  !
  ! count subsets and their degeneracy
  !
  DO ibnd = 2, nbnd
    tdegen = abs( et (ibnd, ik0) - et (ibnd-1, ik0)  ) .lt. eps 
    IF (tdegen) then
      ndeg (nset) = ndeg(nset) + 1
    ELSE
      nset = nset + 1
    ENDIF
  ENDDO
  IF (iverbosity.eq.1) write (stdout, *) & 
     ik0, nset, (ndeg (iset) ,iset=1,nset)
  !
  ! -----------------------------------------------------------
  ! determine rotations within each subset 
  ! and copy into the global matrix
  ! -----------------------------------------------------------
  !

  !$$ initialize the perturbation
  deltav = 0.0d0

  n0 = 0
  DO iset = 1, nset
    !
    !
    !  the size of the small rotation
    nsmall = ndeg (iset)
    !
    ! unimat is initialized to the identity matrix,
    ! so when the subset has dimension 1, just do nothing
    !
    IF (nsmall.gt.1) then
      !
      ! form the matrix elements of the "perturbation"
      !
      allocate ( u (nsmall, nsmall) )

      !$$ allocate matrices for rotation at first
      allocate ( ifail( nsmall ), iwork( 5*nsmall ), w( nsmall ), rwork( 7*nsmall ),&
                 up( nsmall*(nsmall+1)/2 ), cwork( 2*nsmall ), cz( nsmall, nsmall) )

      tdegen = .true.

      !$$ initial seed for perturbation
      na1 = 1
      na2 = 2
      na3 = 3

      na4 = nint(0.5/(omega/dffts%nnr)**(1.0/3.0))
      ! na4 is the number of real space grid points
      ! that roughly corresponds to 0.5 bohr.
      ! (omega is the volume of the unit cell in bohr^3)

      na5 = dfftp%nr1x*dfftp%nr2x/(na4*na4) + dfftp%nr1x/na4
      ! an integer between 0 and na5-1 is going to be picked.
      !$$


      !$$
      !$$ repeat until we find a good u matrix
      !$$
      DO while(tdegen)

         !$$ set up deltav(ir)
         DO i=1,dfftp%nr1x,na4
            DO j=1,dfftp%nr2x,na4
               DO k=1,dfftp%nr3x,na4
                  ir = i + (j - 1) * dfftp%nr1x + (k - 1)*dfftp%nr1x*dfftp%nr2x
                  
                  deltav(ir) = mod(na1*i*j+na2*j*k+na3*k*i,na5)*na4*na4*na4/na5
                  ! here, deltav(ir) is assigned a number so that the
                  ! eigenvalue of the fake perturbation is between 0 and 1, roughly
                  
               ENDDO
            ENDDO
         ENDDO
         !$$

         u = czero
         DO ibnd =1, nsmall
            !
            ! ibnd and jbnd are the indexes of the bands within the subsets
            ! ibndg and jbndg are the global indexes of the bands
            !
            ibndg = ibnd + n0
            !
            aux1(:) = (0.d0, 0.d0)
            DO ig = 1, nng
               aux1 (nls (igktmp (ig) ) ) = evtmp (ig, ibndg)
            ENDDO
            CALL invfft ('Wave', aux1, dffts)
            DO ir = 1, dffts%nnr
               aux1 (ir) = aux1 (ir) * deltav (ir) 
            ENDDO
            CALL fwfft ('Wave', aux1, dffts)
            deltavpsi (1:nng) = aux1 (nls (igktmp (1:nng) ) )
            ! 
            DO jbnd = 1, nsmall
               !
               jbndg = jbnd + n0
               ! 
               u (ibnd, jbnd) = ZDOTC (nng, deltavpsi, 1, evtmp(:, jbndg), 1)
            ENDDO
            !
         ENDDO
         !
         ! ok I veryfied that when deltav(ir)=1, u is the unity matrix (othonormality)
         !
         CALL mp_sum(u, intra_pool_comm)
         !
         ! check hermiticity
         ! 
         DO ibnd = 1, nsmall
            DO jbnd = 1, nsmall
               IF ( abs( conjg (u (jbnd, ibnd)) - u (ibnd, jbnd) ) .gt. eps ) then
                  DO mbnd = 1,nsmall
                     WRITE(stdout,'(10f15.10)') (u (pbnd, mbnd), pbnd=1,nsmall)
                  ENDDO
                  CALL errore ('setphases','perturbation matrix non hermitian',1)
               ENDIF
            ENDDO
         ENDDO
         !
         ! now diagonalize the "perturbation" matrix within the degenerate subset
         !
!$$ allocation is done outside the loop
!$$      allocate ( ifail( nsmall ), iwork( 5*nsmall ), w( nsmall ), rwork( 7*nsmall ),&
!$$                 up( nsmall*(nsmall+1)/2 ), cwork( 2*nsmall ), cz( nsmall, nsmall) )
!$$
         !
         !  packed upper triangular part for zhpevx
         DO jbnd = 1, nsmall
            DO ibnd = 1, nsmall
               up (ibnd + (jbnd - 1) * jbnd/2 ) = u ( ibnd, jbnd)
            ENDDO
         ENDDO
         !
         CALL zhpevx ('V', 'A', 'U', nsmall, up , 0.0, 0.0, &
              0, 0,-1.0, neig, w, cz, nsmall, cwork, &
              rwork, iwork, ifail, info)
         IF (iverbosity.eq.1) then
!$$         if (.true.) then
            !
            WRITE(stdout, '(5x, "Eigenvalues of fake perturbation: ", 10f9.5)') &
                 (w(ibnd),ibnd=1,nsmall)
            WRITE(stdout, *) 
            WRITE(stdout,*) 'before diagonalizing the perturbation:'
            DO ibnd = 1, nsmall
               WRITE(stdout,'(10f9.4)') (u (ibnd, jbnd), jbnd=1, nsmall)
            ENDDO
            !       
            ! recalculate u via similarity transform
            !
            u = matmul ( conjg(transpose( cz)), matmul (u, cz) ) 
            WRITE(stdout,*) 'after diagonalizing the perturbation: via matmul'
            DO ibnd = 1, nsmall
               WRITE(stdout,'(10f9.4)') (u (ibnd, jbnd), jbnd=1, nsmall)
            ENDDO
            WRITE(stdout, '("-----------------"/)') 
            !
         ENDIF
         !
         ! now make sure that all the eigenvalues are nondegenerate
         !
         tdegen = .false.
         DO ibnd = 2, nsmall  
            tdegen = tdegen .or. ( abs( w(ibnd) - w(ibnd-1) ) .lt. eps )  
         ENDDO
         IF(tdegen) write(stdout,*) 'eigenvalues of pert matrix degenerate'

!$$ in the following, instead of killing the program when the perturbation
!$$ gives degenerate eigenvalues, it changes the perturbation and repeat
!$$ until the degeneracy is broken. The perturbation should not be the same
!$$ for wavefunctions with different wavevector or band indices.
!$$
!$$      if (tdegen) call errore ('setphases', &
!$$         'eigenvalues of pert matrix degenerate',1) 
         !
         ! ...that they are nonvanishing... 
         !
!$$         do ibnd = 1, nsmall  
!$$            tdegen = tdegen .or. ( abs( w(ibnd) ) .lt. eps )  
!$$         enddo
!$$      if (tdegen) call errore ('setphases', &
!$$         'eigenvalues of pert matrix too small',1)
         !
         ! ...and that they are not too close to 1 (orthonormality...)
         !
!$$         do ibnd = 1, nsmall  
!$$            tdegen = tdegen .or. ( abs( w(ibnd) - 1.d0 ) .lt. eps )  
!$$         enddo
!$$      if (tdegen) call errore ('setphases', &
!$$         'eigenvalues of pert matrix too close to 1',1) 

         !$$ change the perturbation
         na1 = na1 + 3
         na2 = na2 + 1
         na3 = na3 + 1
         na4 = na4 + 1
         na5 = na5 + 10
         !$$

      ENDDO !$$ repeat until finding a perturbation that gives non-degenerate eigenvalues
      !
      ! copy the rotation for the subset into the global rotation
      !
      n1 = n0 + 1
      n2 = n0 + ndeg(iset)
      unimat ( n1:n2, n1:n2) = cz ( 1:nsmall, 1:nsmall)
      !
      deallocate ( ifail, iwork, w, rwork, up, cwork, cz )
      deallocate ( u )
      !
    ENDIF
    !
    ! update global band index
    !
    n0 = n0 + ndeg(iset)
  ENDDO
  !
  IF (iverbosity.eq.1) then
!$$  if (.true.) then
     WRITE(stdout, *) '--- rotations for unitary gauge ----'
     WRITE(stdout,'(i4)') ik0
     WRITE(stdout,'(8f9.5)') (et(ibnd,ik0),ibnd=1,nbnd)
     DO ibnd = 1, nbnd
       WRITE(stdout,'(10f9.4)') (unimat (ibnd, jbnd), jbnd=1, nbnd)
     ENDDO
  ENDIF
  !
  ! -----------------------------------------------------------
  ! now fix the phases and update rotation matrix
  ! -----------------------------------------------------------
  !
  !  rotate the coefficients with the matrix u 
  ! (to set the phase on the evc's *after* the rotation)
  !
  cnew(:,:) = czero
  DO ibnd = 1, nbnd
    DO jbnd = 1, nbnd
      cnew(1:nnglen, ibnd) = cnew(1:nnglen, ibnd) & 
        + unimat (jbnd, ibnd) * c (1:nnglen, jbnd)
    ENDDO
  ENDDO
  !
  ! for every band, find the largest coefficient and determine
  ! the rotation which makes it real and positive 
  !
  DO ibnd = 1, nbnd
    ! 
    !  this is to identify the largest c(G) by using *integer* 
    !  comparisons on ndig digits [see remarks at the bottom about this]
    !
    imaxnorm = 0
    DO ig = 1, nnglen
      ctmp = cnew (ig, ibnd)
      tmp = REAL(conjg ( ctmp ) * ctmp)
      itmp = nint (10.d0**dble(ndig) * tmp)
      IF (itmp.gt.imaxnorm) then
        imaxnorm = itmp
        igmax = ig
      ENDIF
    ENDDO
    !
    ctmp = cnew (igmax, ibnd)
    tmp = REAL(conjg ( ctmp ) * ctmp)
    !
    ! ...and the corresponding phase
    !
    IF ( abs(tmp) .gt. epsc ) then
      ! note that if x + i * y = rho * cos(theta) + i * rho * sin(theta),
      ! then theta = atan2 ( y, x) (reversed order of x and y!)
      theta = atan2 ( aimag( ctmp ), real ( ctmp ) )
    ELSE
      CALL errore ('setphases','cnew = 0 for some bands: increase nnglen',1)
    ENDIF
    !  
    IF (iverbosity.eq.1) then
!$$    if (.true.) then
      WRITE(stdout, '(3i4,2x,f15.10,2(3x,2f9.5))') ik0, ibnd, igmax, theta/pi*180.d0, &
         ctmp, ctmp * exp ( -ci * theta), exp ( -ci * theta)
    ENDIF
    !
    !  now cancel this phase in the rotation matrix
    !
    unimat (:, ibnd) = unimat (:, ibnd) * exp ( -ci * theta)
    !
  ENDDO
  !
  IF (iverbosity.eq.1) then
!$$  if (.true.) then
     WRITE(stdout, *) '--- rotations including phases  ----'
     DO ibnd = 1, nbnd
       WRITE(stdout,'(10f9.4)') (unimat (ibnd, jbnd), jbnd=1, nbnd)
     ENDDO
  ENDIF
  !
  !  last check: unitarity
  !  (if this test is passed, even with wrong phases the final
  !  results -nonwannierized- should be ok)
  !
  unitcheck = matmul ( conjg( transpose (unimat)), unimat)
  DO ibnd = 1, nbnd
    DO jbnd = 1, nbnd
      IF ( (ibnd.ne.jbnd) .and. ( abs(unitcheck (ibnd, jbnd)) .gt. eps )                  &
     .or.  (ibnd.eq.jbnd) .and. ( abs ( abs(unitcheck (ibnd, jbnd)) - 1.d0 ) .gt. eps ) )&
        CALL errore ('setphases','final transform not unitary',1)
    ENDDO
  ENDDO
  !
  nullify ( evtmp ) 
  nullify ( igktmp ) 
  !
  end subroutine setphases
!---------------------------------------------------------------------
!
! NOTA BENE: I truncate the c(G) to a given number of digits in
! order to guarantee that the same norm-ordering of c(G) holds
! for different q-runs, machines, libraries etc:
! run q  = 0 : ig = 3 |c(igmax)| = 0.34.....263  <-- igmax = 3
! run q  = 0 : ig = 4 |c(igmax)| = 0.34.....261
! run q /= 0 : ig = 3 |c(igmax)| = 0.34.....260
! run q /= 0 : ig = 4 |c(igmax)| = 0.34.....265  <-- igmax = 4
! In the situation above, I will have psi(r) with a phase difference
! corresponding to that of c(3) and c(4)...
! (Note that the ordering of the G-vectors is always the same) 
! Mind : ndig should be smaller than machine precision, otherwise the 
! integere comparison is useless
!
!---------------------------------------------------------------------
!
