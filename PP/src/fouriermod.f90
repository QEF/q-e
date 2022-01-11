!
! Copyright (C) 2001-2022 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
! Author: Ivan Carnimeo (September 2021)
!----------------------------------------------------------------------------
!
!----------------------------------------------------------------------------
MODULE fouriermod
  USE kinds,                ONLY : dp
implicit none
save
  real(dp), parameter :: eps = 0.000010d0, Zero = 0.0d0, One = 1.0d0, Two = 2.0d0, Four = 4.0d0
  real(dp), parameter :: Pi = Four*atan(One)
  !
  ! the largest Miller index used to generate all the lattice vectors inside an outer shell
  integer :: NMax
  integer :: NC  
  real(dp), allocatable :: C(:) 
  !
  integer :: NStars                       ! total number of Star functions 
  real(dp), allocatable :: VecStars(:,:)  ! symmetry inequivalent lattice vectors generating the Star functions (one per Star)
  integer :: NUser                        ! (Optional) number of user-given star functions 
  real(dp), allocatable :: VecUser(:,:)   ! (Optional) user-given star functions 
  !
CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine allocate_fourier( )
implicit none
  allocate ( C(NC) ) 
  return
end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine fourier(Nb, Nq, q, eq, Nk, k, ek, Nsym, at, bg, Op)
!
! compute the band structure with Fourier interpolation from
! D. D. Koelling, J. H. Wood, J. Comput. Phys., 67, 253-262 (1986)
!
implicit none
  integer,  intent(in) :: Nq, Nk, Nb, Nsym
  real(dp), intent(in) :: q(3,Nq), k(3,Nk), eq(Nq,Nb)
  real(dp), intent(in) :: at(3,3) 
  real(dp), intent(in) :: bg(3,3) 
  real(dp), intent(in) :: Op(1:3,1:3,1:Nsym)
  real(dp), intent(out) :: ek(Nk,Nb)
  !
  ! local variables
  !
  complex(dp), allocatable :: fStarsOnQ(:,:) ! Star functions at uniform q-points
  complex(dp), allocatable :: fStarsOnK(:,:) ! Star functions at path of k-points 
  complex(dp), allocatable :: matQQ(:,:)     ! this is exactly H in the reference article
  complex(dp), allocatable :: matKQ(:,:)     ! this is an intermediate quantity S_m(q)*S_m(k)/rho(R_m) to construct J
  complex(dp), allocatable :: matJ(:,:)      ! this is exactly J in the reference article
  complex(dp), allocatable :: ek_c(:,:), eq_c(:,:) ! complex band energies (for ZGEMM)
  !
  ! for matrix inversion
  complex(dp), allocatable :: matQQ_(:,:)     
  real(dp), allocatable :: rmatQQ(:,:)     
  real(dp), allocatable :: rmatQQ_(:,:)     
  real(dp), allocatable :: rmatJ(:,:)     
  !
  write(*,'(A)') '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
  write(*,'(A)') 'Fourier interpolation method'
  write(*,*) 'Creating Star functions...'
  Call find_stars(NSym, Op, at) 
  !
  write(*,*) 'Checking Star functions periodicity...'
  Call check_stars(Nq, q, NSym, Op, bg) 
  !
  ! fStarsOnQ = S_m(q) / sqrt(rho_m)
  write(*,*) 'Computing fStarsOnQ...'
  allocate( fStarsOnQ(Nq,NStars) )
  fStarsOnQ = (Zero, Zero)
  Call compute_stars(fStarsOnQ, Nq, Nq, q, NSym, Op, 2) 
  !
!!!!civn 
!!!  ! matQQ = fStarsOnQ * fStarsOnQ^T  = sum_i S_m(q_i) S_n(q_i) 
!!!  write(*,*) 'Checking ortogonality fStarsOnQ^T * fStarsOnQ...'
!!!  allocate(matQQ(NStars,NStars), rmatQQ(NStars,NStars))
!!!  matQQ = (Zero, Zero)
!!!  Call ZGEMM('C', 'N', NStars, NStars, Nq, (One,Zero), fStarsOnQ, Nq, fStarsOnQ, Nq, (Zero,Zero), matQQ, NStars)
!!!  matQQ = matQQ / dble(Nq)
!!!  Call matprt_k('matQQ', NStars, NStars, matQQ)
!!!  rmatQQ = dble(matQQ)
!!!  Call MatCheck('rmatQQ',rmatQQ,NStars,NStars)
!!!  deallocate(matQQ, rmatQQ)
!!!stop  
  !
  ! matQQ = fStarsOnQ * fStarsOnQ^T  = sum_m S_m(q_i) S_m(q_j) / rho_m
  write(*,*) 'Computing fStarsOnQ * fStarsOnQ^T...'
  allocate(matQQ(Nq,Nq))
  matQQ = (Zero, Zero)
  Call ZGEMM('N', 'C', Nq, Nq, NStars, (One,Zero), fStarsOnQ, Nq, fStarsOnQ, Nq, (Zero,Zero), matQQ, Nq)
  !
  ! matQQ --> matQQ^(-1) 
  write(*,*) 'Inverting fStarsOnQ * fStarsOnQ^T...'
  allocate( rmatQQ(Nq,Nq), rmatQQ_(Nq,Nq), rmatJ(Nq,Nq) )
  rmatQQ(:,:)  = dble(matQQ(:,:))
  rmatQQ_(:,:) = rmatQQ(:,:) 
  rmatJ(:,:)   = Zero 
  Call MatInv('G', Nq, rmatQQ)
  write(*,*) 'Checking inverse...'
  Call DGEMM("N", "N", Nq, Nq, Nq, One, rmatQQ, Nq, rmatQQ_, Nq, Zero, rmatJ, Nq)
  Call MatCheck('rmatJ',rmatJ,Nq,Nq)
  matQQ(:,:) = (One,Zero) * rmatQQ(:,:)
  deallocate( rmatQQ, rmatQQ_, rmatJ )
  !
  ! fStarsOnK = S_m(k) / sqrt(rho_m)
  write(*,*) 'Computing fStarsOnK...'
  allocate( fStarsOnK(Nk,NStars) )
  fStarsOnK = (Zero, Zero)
  Call compute_stars(fStarsOnK, Nk, Nk, k, NSym, Op, 2) 
  !
  ! matKQ = fStarsOnK * fStarsOnQ^T
  allocate(matKQ(Nk,Nq))
  matKQ = (Zero, Zero)
  Call ZGEMM('N', 'C', Nk, Nq, NStars, (One,Zero), fStarsOnK, Nk, fStarsOnQ, Nq, (Zero,Zero), matKQ, Nk)
  !
  ! matJ = matKQ * matQQ^(-1)
  allocate(matJ(Nk,Nq))
  matJ = (Zero, Zero)
  Call ZGEMM('N', 'N', Nk, Nq, Nq, (One,Zero), matKQ, Nk, matQQ, Nq, (Zero,Zero), matJ, Nk)
  !
  ! e(k) = matJ * e(q)
  allocate(ek_c(Nk,Nb), eq_c(Nq,Nb))
  eq_c(:,:) = (One,Zero) * eq(:,:)
  ek_c(:,:) = (Zero,Zero)
  Call ZGEMM('N', 'N', Nk, Nb, Nq, (One,Zero), matJ, Nk, eq_c, Nq, (Zero,Zero), ek_c, Nk)
  !
  ek(:,:) = dble(ek_c(:,:))
  !
  deallocate( C ) 
  deallocate( ek_c, eq_c )
  deallocate( VecStars )
  deallocate( matQQ )
  deallocate( matKQ )
  deallocate( matJ )
  deallocate( fStarsOnQ )
  deallocate( fStarsOnK  )
  !
  return
  !
end subroutine fourier 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine find_stars(NSym, Op, at, Skip000_)
implicit none
  integer,  intent(in) :: NSym
  real(dp), intent(in) :: at(3,3)  
  real(dp), intent(in) :: Op(1:3,1:3,1:Nsym)
  logical, optional, intent(in) :: Skip000_ ! whether to skip Miller indeces (0 0 0)
  !
  ! local variables
  logical :: Skip000 
  !
  ! arrays containing the lattice vectors inside the outer shell defined by NMax
  integer :: NAll, NPrint
  real(dp), allocatable :: VecAll(:,:) ! all lattice vectors
  real(dp), allocatable :: VecInq(:,:) ! symmetry inequivalent lattice vectors
  real(dp), allocatable :: ModAll(:)   ! modules of the lattice vectors
  integer, allocatable  :: MapAll(:)   ! permutation map of the lattice vectors
  !
  integer :: ii, jj, kk, isym, ivec, jvec
  real(dp) :: rdiff
  real(dp) :: vec_i(3), mod_i, vecOp_i(3)
  real(dp) :: vec_j(3), mod_j
  !
  if(present(Skip000_)) then 
    Skip000 = Skip000_
  else
    Skip000 = .false.
  endif 
  !
  NAll = (2 * NMax + 1 )**3 ! from -Nmax to Nmax is (2*Nmax + 1), for 3 space directions 
  if(Skip000) then 
    NAll = NAll - 1  ! remove the (0, 0, 0) lattice vector 
    write(*,*) 'Skipping the (0,0,0) lattice vector...'
  else
    write(*,*) 'Including the (0,0,0) lattice vector...'
  end if 
  !
  if(NUser.gt.0) then 
    NAll = NAll + NUser  
    write(*,*) "Creating ", NAll, " vectors from ", NMax, " indexes and ", NUser, " user-given vectors"
  else
    write(*,*) "Creating ", NAll, " vectors from ", NMax, " indexes"
  end if
  !
  allocate ( VecAll(3,NAll), ModAll(NAll), MapAll(NAll) ) 
  !
  ivec = 0 
  do ii = -NMax, NMax
    do jj= -NMax, NMax
      do kk= -NMax, NMax
        if(Skip000.and.(ii.eq.0.and.jj.eq.0.and.kk.eq.0)) cycle 
        ivec = ivec + 1 
        VecAll(:,ivec) = dble(ii)*at(:,1) + dble(jj)*at(:,2) + dble(kk)*at(:,3)   
        ModAll(ivec) = sqrt( VecAll(1,ivec)**2 + VecAll(2,ivec)**2 + VecAll(3,ivec)**2 ) 
        MapAll(ivec) = ivec
      end do 
    end do 
  end do 
  !
  if(NUser.gt.0) then 
    do ii = 1, NUser
      ivec = ivec + 1  
      VecAll(:,ivec) = VecUser(:,ii)
      ModAll(ivec) = sqrt( VecAll(1,ivec)**2 + VecAll(2,ivec)**2 + VecAll(3,ivec)**2 ) 
      MapAll(ivec) = ivec
    end do 
  end if 
  !
  if(ivec.ne.NAll) then 
    write(*,*) "ERROR: wrong number of lattice vectors for a given NMax"
    write(*,*) "NMax= ",NMax," ivec=",ivec," NAll=",NAll, " NUser=",NUser
    stop
  endif
  !  
  write(*,*) "Sorting vectors in shells..."
  Call hpsort_eps (NAll, ModAll, MapAll, eps)
  !
  write(*,*) "Removing symmetry equivalent lattice vectors..." 
  if(NUser.gt.0) write(*,*) "WARNING: user-given vectors will be removed they are symmetry equivalent"
  !
  allocate( VecInq(3,NAll) ) 
  VecInq = Zero 
  !
  NPrint = NAll / 10
  NStars = 0 
  do jj = 1, NAll
    !
    if(jj.eq.1.or.mod(jj,NPrint).eq.0) write(*,'(5X,I10,A,f12.2,A,I10,A)') &
              jj-1, ' (',dble(100*(jj-1))/dble(NAll), '%) vectors analysed ... ', NStars, " Stars found"
    !
    jvec = MapAll(jj)
    vec_j(:) = VecAll(:,jvec) ! this is in the original ordering (jvec)
    mod_j = ModAll(jj) ! sqrt(dot_product(vec_j,vec_j)) ! this has been sorted by hpsort_eps (jj)
    !
    ! this loop is needed because in principle there could be two Stars
    ! generated by different vectors with the same module 
    do ivec = 1, NStars 
      vec_i(:) = VecInq(:,ivec)
      mod_i = sqrt(dot_product(vec_i,vec_i))
      !
      if(abs(mod_i-mod_j).lt.eps) then ! if the modules are different they cannot be equivalent
        do isym = 1, NSym
          Call applyOp(isym, Op(1:3,1:3,isym), vec_i, vecOp_i)
          rdiff=sqrt( (vecOp_i(1)-vec_j(1))**2 + (vecOp_i(2)-vec_j(2))**2 + (vecOp_i(3)-vec_j(3))**2 ) 
          if(rdiff.lt.eps) go to 10 ! vec_j is equivalent to vec_j, skip it
        end do 
      end if 
      !
    end do  
    !
    ! vec_j is not equivalent to any vec_i : it generates a new Star
    NStars = NStars + 1
    VecInq(:,NStars) = vec_j(:)
    !
10  continue
    !
  end do 
  !
  write(*,*) NStars, " Stars of symmetry inequivalent lattice vectors found..."
  !
  ! Knowing NStars we can now allocate VecStars with the right size 
  ! and deallocate the over-sized VecInq
  !
  allocate( VecStars(3,NStars) )
  !
  VecStars = Zero
  !
  VecStars(:,:) = VecInq(:,1:NStars)
  !
  deallocate( VecAll, ModAll, MapAll, VecInq )
  if(allocated(VecUser)) deallocate(VecUser)
  !
  return 
  !
end subroutine find_stars
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine check_stars(Np, p, NSym, Op, bg)
! Check whether the Star functions have the periodicity 
! of the reciprocal lattice
implicit none
  integer, intent(in) :: Np, NSym
  real(dp), intent(in) :: p(1:3,1:Np), Op(1:3,1:3,1:NSym)
  real(dp), intent(in) :: bg(3,3)  
  real(dp), parameter  :: Thr=0.0010d0
  real(dp) :: vec(3), gvec(3), pvec(3)
  integer :: istar, ip, ig, jg, kg, isym
  complex(dp) :: fp, fpg 
  !
  do istar = 1, NStars
    vec(:) = VecStars(:,istar)
    do ip = 1, Np 
      fp = star_function(0, p(1:3,ip), vec, NSym, Op)
      do ig = -1, 1
        do jg = -1, 1
          do kg = -1, 1
            gvec(:) = dble(ig) * bg(:,1) + dble(jg) * bg(:,2) + dble(kg) * bg(:,3)
            pvec(:) = p(:,ip) + gvec(:)
            fpg = star_function(0, pvec(:), vec, NSym, Op)
            if(abs(fp-fpg).gt.Thr) then 
              if(NUser.gt.0) then 
                write(*,'(A)') 'WARNING: A Star function does not have reciprocal lattice periodicity'
              else
                write(*,'(A)') 'ERROR: A Star function does not have reciprocal lattice periodicity'
              end if 
              write(*,'(A,I5,A,3f12.4)')           'istar:   ', istar, '    vec: ', vec(:)
              write(*,'(A,I5,A,3f12.4,A,2f24.12)') 'ip:      ', ip,    '  P-vec: ', p(:,ip), ' fp: ', fp  
              write(*,'(3(A,I5),A,3f12.4)')        'ig: ', ig, ' jg: ', jg, ' kg: ', kg, ' G-vec: ', gvec(:)
              write(*,'(A,3f12.4,A,2f24.12)')      'P+G-vec: ', pvec(:), ' fpg: ', fpg
              if(NUser.gt.0) then 
                go to 30
              else
                stop
              end if
            end if 
          end do 
        end do 
      end do 
      isym = NSym/3 + 1 
      pvec = Zero
      Call applyOp(isym, Op(1:3,1:3,isym), p(1:3,ip), pvec)
      fpg = star_function(0, pvec(:), vec, NSym, Op)
      if(abs(fp-fpg).gt.Thr) then
        if(NUser.gt.0) then 
          write(*,'(A)') 'WARNING: A Star function does not have reciprocal lattice periodicity'
        else
          write(*,'(A)') 'ERROR: A Star function does not have reciprocal lattice periodicity'
        end if 
        write(*,'(A,I5,A,3f12.4)')           'istar:   ', istar, '    vec: ', vec(:)
        write(*,'(A,I5,A,3f12.4,A,2f24.12)') 'ip: ',    ip,   ' P-vec: ', p(:,ip), ' fp: ',  fp  
        write(*,'(A,I5,A,3f12.4,A,2f24.12)') 'isym:  ', isym, ' OpP: ',   pvec(:), ' fpg: ', fpg
        if(NUser.le.0) stop
      end if
    end do 
30  continue   
  end do 
  !
  return
  !
end subroutine check_stars
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine compute_stars(A, LDA, Np, p, NSym, Op, ialpha, DoDiff_, S)
!
! Computes the matrix A(LDA, NStars):
!
!         A_im = alpha * S_m(p_i)
!
! where: 
!         m (=istar) is the index of a Star vector
!         p is a set of Np k-points  (either q or k) 
!         i = 1...LDA and LDA=Np
!         alpha = One               ... ialpha = 0 
!                 1/rho(R_m)        ... ialpha = 1 
!                 1/sqrt(rho(R_m))  ... ialpha = 2 
!
! if DoDiff_ is true, then the matrix is:
!
!         A_im = alpha * [ S_m(p_i) - S_m(p_Np) ] 
!
! where: 
!         i = 1...LDA and LDA=Np-1
! and the array S_m(p_Np) is also returned as output
!
implicit none
  complex(dp), intent(out) :: A(LDA,NStars)
  complex(dp), optional, intent(out) :: S(NStars)
  integer, intent(in) :: LDA, Np, NSym, ialpha
  real(dp), intent(in) :: p(1:3,1:Np), Op(1:3,1:3,1:NSym)
  logical, optional, intent(in) :: DoDiff_
  !
  logical :: DoDiff
  !
  real(dp) :: vec(3)  ! R_m
  real(dp) :: sqrtrho ! sqrt(rho(R_m))
  integer :: istar, ip
  complex(dp) :: alpha 
  complex(dp) :: fStar  ! S_m(p_i)
  complex(dp) :: fStarN ! S_m(p_Np)
  !
  if(ialpha.lt.0.or.ialpha.gt.2) then 
    write(*,*) 'ERROR: wrong ialpha in compute_stars'
    stop
  end if
  !
  if(present(DoDiff_)) then 
    DoDiff = DoDiff_
    if(.not.present(S)) then 
      write(*,*) 'ERROR: please provide S with DoDiff=.true.'
      stop
    endif 
  else
    DoDiff = .false.
  end if 
  !
  write(*,'(A,L,3(A,I5))') 'compute_stars: DoDiff: ', DoDiff, ' LDA: ', LDA, ' Np:', Np, ' NStars:', NStars
  if((DoDiff.and.(LDA.ne.(Np-1))).or.(.not.DoDiff.and.(LDA.ne.Np))) then  
    write(*,*) 'ERROR: Wrong dimensions in compute_stars'
    stop
  end if
  !
  do istar = 1, NStars
    vec(:) = VecStars(:,istar)  ! R_m
    if(ialpha.eq.0) then 
      alpha  = (One,Zero) 
    elseif(ialpha.eq.1) then 
      sqrtrho = sqrt_rho(vec)
      alpha  = (One,Zero) / sqrtrho / sqrtrho
    elseif(ialpha.eq.2) then 
      alpha = (One,Zero) / sqrt_rho(vec)
    end if
    fStarN  = (Zero, Zero)
    if(DoDiff) then 
      fStarN  = star_function(0, p(1:3,Np), vec, NSym, Op)
      S(istar) = fStarN
    end if 
    do ip = 1, LDA
      fStar = star_function(0, p(1:3,ip), vec, NSym, Op) 
      A(ip,istar) = alpha * (fStar-fStarN) 
    end do 
    if(abs(aimag(A(ip,istar))).gt.eps) write(*,*) "Star function: ", ip, istar, A(ip,istar), " WARNING non zero imaginary part!!"
  end do 
  !
  return
  !
end subroutine compute_stars
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
complex(dp) function star_function (iprint, p, vec, NSym, Op)
implicit none
  integer, intent(in) :: iprint
  integer, intent(in) :: NSym
  real(dp), intent(in) :: p(1:3), vec(3)
  real(dp), intent(in) :: Op(1:3,1:3,1:NSym)
  ! local variables
  real(dp) :: vecOp(3)
  real(dp) :: diffMod 
  real(dp), allocatable :: vecInq(:,:) ! vecInq(NSym,3) a list of inequivalent vectors to remove duplicates
  complex(dp) :: carg, cfunc
  integer :: isym, jsym, NVec
  !
  cfunc = (Zero, Zero)
  !
  allocate( vecInq(3,NSym) )
  !
  vecInq(:,:) = Zero 
  NVec = 0 
  !
  do isym = 1, NSym
    Call applyOp(isym, Op(1:3,1:3,isym), vec, vecOp)
    do jsym = 1, NVec
      ! for vecInq=Zero diffMod never .lt.eps
      diffMod = sqrt((vecOp(1)-vecInq(1,jsym))**2 + (vecOp(2)-vecInq(2,jsym))**2 + (vecOp(3)-vecInq(3,jsym))**2 )
      if(iprint.gt.0) write(*,'(2(I5,3f6.2,3x))') isym, vecOp(:), jsym, vecInq(:,jsym)
      if(diffMod.lt.eps) go to 20 ! vecOp already included
    end do 
    ! vecOp is new
    NVec = NVec + 1 
    vecInq(:,NVec) = vecOp(:)
    carg  = (Zero,One) * Two * Pi *  dot_product(vecOp,p) 
    cfunc = cfunc + exp(carg)
20  continue
    if(iprint.gt.0) write(*,'(I5,3x,2f12.4)') NVec, cfunc
  end do 
  carg = (One,Zero) * sqrt(dble(NVec))
  star_function = cfunc/carg 
  !
  deallocate( vecInq )
  !
  return
  !
end function star_function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real(dp) function sqrt_rho(vec)
! computes:
!           sqrt(C0 + C1 * |vec|^2 + C2 * |vec|^4)
implicit none
  real(dp), intent(in) :: vec(1:3) 
  real(dp) :: rmod, rho  
  integer :: ic, iexp
  !
  rho = C(1)
  if(NC.gt.1) then 
    rmod = sqrt(dot_product(vec,vec))
    do ic  = 2, NC
      iexp = 2*(ic - 1 )
      rho = rho + C(ic) * rmod**iexp
    end do 
  end if 
  sqrt_rho = sqrt(rho)
  !
  return  
  !
end function sqrt_rho
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine applyOp(isym, OpMat, vec, vecOp)
! apply symmetry operation in OpMat:
!         vecOp = OpMat * vec
implicit none 
  integer, intent(in) :: isym
  real(dp), intent(in) :: OpMat(1:3,1:3), vec(1:3)
  real(dp), intent(out) :: vecOp(1:3)
  real(dp) :: vecErr
  !
  vecOp(:) = matmul(OpMat, vec)    
  !
  vecErr = abs(sqrt(dot_product(vec,vec))-sqrt(dot_product(vecOp,vecOp)))
  !
  if (vecErr.gt.eps) then   
    write(*,*) "ERROR: non-unitary symmetry operation found"
    write(*,*) "isym: ", isym
    write(*,*) "vec:     ", vec(:)
    write(*,*) "vecOp:   ", vecOp(:)
    write(*,*) "vecErr:  ", vecErr
    Call MatPrt("OpMat", 3, 3, OpMat)
    stop
  endif
  !
  return
  !
end subroutine applyOp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE
