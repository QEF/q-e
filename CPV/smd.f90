!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!*******************************************************************
!
!       String Method subroutines for implementation with
!       O-sesame/CP code
!
!       Y. Kanai
!
!       DATE    : 050504
!
!
! sub. sminit
! sub. tangent
! sub. perp
! sub. linearp
! sub. arc
! sub. absvec
! sub. startcoord
! sub. rot
! sub. tran
! sub. distatoms
! sub. init_path
!
!
!
!-----------------------------------------------------------------------
subroutine sminit (ibrav,celldm, ecut, ecutw,ndr,nbeg,  &
     tfirst,delt,tps, iforce)
  !-----------------------------------------------------------------------
  !
  !     initialize G-vectors and related quantities
  !     use ibrav=0 for generic cell vectors given by the matrix h(3,3)
  !
  use control_flags, only: iprint, thdyn
  use io_global, only: stdout
  use mp_global, only: nproc_image
  !use gvec
  use gvecw, only: ngw
  use ions_base, only: na, pmass, nsp, randpos, nat
  use cell_base, only: ainv, a1, a2, a3, r_to_s, s_to_r
  use electrons_base, only: nx => nbspx, f, nudx, nspin
  use constants, only: pi, fpi
  use cell_base, only: hold, h
  use betax, only: mmx, refg
  !use restartsm
  use cp_interfaces, only: writefile, readfile
  use parameters, only: nacx, nhclm
  USE smd_rep, only: rep
  USE path_variables, only: &
        sm_p => smd_p
  USE fft_base,         ONLY: dfftb, fft_dlay_descriptor
  USE fft_types,        ONLY: fft_box_allocate


  implicit none
  ! input/output
  integer ibrav, ndr, nbeg, iforce(3,nat)
  logical tfirst
  real(8) celldm(6), ecut, ecutw
  real(8) delt
  ! local
  real(8) randy
  integer i, j, ia, is, nfi, isa, isat


  ! YK
  ! present in the call to read(p)file, not actually used
  !      complex(8) c0(1,1),cm(1,1)
  !      real(8) taum(1,1,1),vel(1,1,1),velm(1,1,1),acc(nacx)
  !      real(8) lambda(1,1),lambdam(1,1)
  !      real(8) xnhe0,xnhem,vnhe,xnhp0,xnhpm,vnhp, ekincm
  !      real(8) xnhh0(3,3),xnhhm(3,3),vnhh(3,3),velh(3,3)
  !      real(8) fion(1,1,1)

  ! YK
  complex(8) c0(ngw,nx),cm(ngw,nx)
  real(8) taum(3,nat),vel(3,nat),velm(3,nat),acc(nacx)
  real(8) lambda(nudx,nudx,nspin),lambdam(nudx,nudx,nspin)
  real(8) xnhe0,xnhem,vnhe,xnhp0(nhclm),xnhpm(nhclm),vnhp(nhclm), ekincm
  real(8) xnhh0(3,3),xnhhm(3,3),vnhh(3,3),velh(3,3)
  real(8) fion(3,nat),tps
  real(8) mat_z(1,1,1)
  integer      nhpcl, nhpdim
  !


  integer :: sm_k,smpm 
  !
  !
  smpm = sm_p -1
  !
  ! taus = scaled, tau0 = alat units
  !
  DO sm_k=0, sm_p
     call r_to_s( rep(sm_k)%tau0, rep(sm_k)%taus, na, nsp, ainv )
  ENDDO
  !
  !
  !  Allocate box descriptor
  !
  CALL fft_box_allocate( dfftb, nproc_image, nat )
  !
  !
  WRITE( stdout,*) '   NOTA BENE: refg, mmx = ', refg, mmx
  !
  if( nbeg >= 0 ) then
     !
     ! read only h and hold from file ndr
     !
     call readfile                                              &
          &     (-1,ndr+1,h,hold,nfi,c0,cm,rep(0)%tau0,taum,vel,velm,acc,   &
          &       lambda,lambdam,xnhe0,xnhem,vnhe,xnhp0,xnhpm,vnhp,nhpcl,nhpdim,ekincm,   &
          &       xnhh0,xnhhm,vnhh,velh,ecut,ecutw,delt,pmass,ibrav,celldm, &
          &       fion,tps, mat_z, f)
     !
     WRITE( stdout,344) ibrav
     do i=1,3
        WRITE( stdout,345) (h(i,j),j=1,3)
     enddo
     WRITE( stdout,*)

  else
     !
     ! with variable-cell we use h to describe the cell
     !
     do i = 1, 3
        h(i,1) = a1(i)
        h(i,2) = a2(i)
        h(i,3) = a3(i)
     enddo

     hold = h

  end if
  !
  !     ==============================================================
  !     ==== generate true g-space                                ==== 
  !     ==============================================================
  !
  call newinit( h )
  !
344 format(' ibrav = ',i4,'       cell parameters ',/)
345 format(3(4x,f10.5))

  return
end subroutine sminit

! ======================================================================!

!
!================================================================
!
!       subroutine TANGENT
!
!       TANGENT subroutine for O-sesame/CP implementation
!
!
!       REF. Upwind Scheme of JCP 113, 9978 (2000) 
!
!================================================================

SUBROUTINE TANGENT(state,tan)


  use ions_base, ONLY: na, nsp, nat
  use smd_ene, ONLY: etot_ar
  USE path_variables, only: &
        sm_p => smd_p, &
        ptr => smd_ptr

  IMPLICIT NONE

  integer :: i,j,is,ia,sm_k,smpm, isa

  type(ptr) :: state(0:sm_p)
  type(ptr) :: tan(0:sm_p)

  real(8) :: ene1,ene2,tmp1,tmp2
  real(8), parameter :: zero = 1.d-15

  ! ---------------------------- !


  tan(0)%d3(:,1:nat) = 0.d0
  tan(sm_p)%d3(:,1:nat) = 0.d0   

  smpm = sm_p -1


  DO sm_k=1,smpm      ! >>>>>>>>>>>>>>>>>>>>>>  !


     IF((etot_ar(sm_k+1) >= etot_ar(sm_k)) .AND. &
          &  (etot_ar(sm_k) >= etot_ar(sm_k-1))) THEN

        isa = 0
        DO is=1,nsp
           DO ia=1,na(is)
              isa = isa + 1
              DO i=1,3
                 tan(sm_k)%d3(i,isa)=state(sm_k+1)%d3(i,isa) & 
                      & -state(sm_k)%d3(i,isa)
              ENDDO
           ENDDO
        ENDDO

     ELSE IF((etot_ar(sm_k+1) <= etot_ar(sm_k)) .AND. &
          & (etot_ar(sm_k) <= etot_ar(sm_k-1))) THEN

        isa = 0
        DO is=1,nsp
           DO ia=1,na(is)
              isa = isa + 1
              DO i=1,3
                 tan(sm_k)%d3(i,isa)=state(sm_k)%d3(i,isa) &
                      & -state(sm_k-1)%d3(i,isa)
              ENDDO
           ENDDO
        ENDDO

     ELSE

        tmp1=ABS(etot_ar(sm_k-1)-etot_ar(sm_k))
        tmp2=ABS(etot_ar(sm_k)-etot_ar(sm_k+1))
        ene1=MAX(tmp1,tmp2)
        ene2=MIN(tmp1,tmp2)

        IF(etot_ar(sm_k+1) >= etot_ar(sm_k-1)) THEN

           isa = 0
           DO is=1,nsp
              DO ia=1,na(is)
                 isa = isa + 1
                 DO i=1,3

                    tan(sm_k)%d3(i,isa) = ene1*(state(sm_k+1)%d3(i,isa) & 
                         & -state(sm_k)%d3(i,isa)) &
                         & + ene2*(state(sm_k)%d3(i,isa)  &
                         & -state(sm_k-1)%d3(i,isa))

                 ENDDO
              ENDDO
           ENDDO

        ELSE

           isa = 0
           DO is=1,nsp
              DO ia=1,na(is)
                 isa = isa + 1
                 DO i=1,3

                    tan(sm_k)%d3(i,isa) = ene2*(state(sm_k+1)%d3(i,isa) &
                         & -state(sm_k)%d3(i,isa)) &
                         & + ene1*(state(sm_k)%d3(i,isa)  &
                         & -state(sm_k-1)%d3(i,isa))

                 ENDDO
              ENDDO
           ENDDO

        ENDIF

     ENDIF


     ! Normalization of tangent ------------------------------- 

     tmp1=0.d0

     isa = 0
     DO is=1,nsp
        DO ia=1,na(is)
           isa = isa + 1
           DO i=1,3

              tmp1 = tmp1 + tan(sm_k)%d3(i,isa)**2.d0

           ENDDO
        ENDDO
     ENDDO


     tmp1 = SQRT(tmp1)

     isa = 0
     DO is=1,nsp
        DO ia=1,na(is)
           isa = isa + 1
           DO i=1,3

              IF(ABS(tan(sm_k)%d3(i,isa)) >= zero ) THEN
                 tan(sm_k)%d3(i,isa) = tan(sm_k)%d3(i,isa)/tmp1
              ELSE
                 tan(sm_k)%d3(i,isa)=0.d0
              ENDIF

           ENDDO
        ENDDO
     ENDDO


  ENDDO           !  <<<<<<<<<<<<< SM LOOP <<<<<<    ! 

  RETURN


END SUBROUTINE TANGENT


!=======================================================================
!
! subroutine PERP       for O-sesame/CP implementation
!
!
!       This calculates the vector perpendicular to the tangent for
!       a given vector.
!
!======================================================================

SUBROUTINE PERP(vec,tan,paraforce)


  use ions_base, ONLY: na, nsp, nat

  IMPLICIT NONE

  integer :: i,is,ia,isa

  real(8) :: vec(3,nat),tan(3,nat)
  real(8) :: paraforce
  real(8) :: dotp

  !----------------------------------------------


  ! calculate DOT product of the vectors -----

  dotp = 0.d0

  isa = 0
  DO is=1,nsp
     DO ia=1,na(is)
        isa = isa + 1
        DO i=1,3
           dotp = dotp + vec(i,isa)*tan(i,isa)
        ENDDO
     ENDDO
  ENDDO

  paraforce = dotp


  ! calculate PERP of the vec  -----------

  isa = 0
  DO is=1,nsp
     DO ia=1,na(is)
        isa = isa + 1
        DO i=1,3
           vec(i,isa) = vec(i,isa) - dotp*tan(i,isa)
        ENDDO
     ENDDO
  ENDDO

  RETURN

END SUBROUTINE PERP


!=======================================================================
!
! subroutine LINEARP 
!
!  calculates the linear interpolated replicas between
!  two replicas.
!
!
!=======================================================================

SUBROUTINE LINEARP(state)


  use ions_base, ONLY: na, nsp
  USE path_variables, only: &
        sm_p => smd_p, &
        ptr => smd_ptr


  IMPLICIT NONE

  integer :: i,j,is,ia,sm_k,smpm,isa

  type(ptr) :: state(0:sm_p)

  real(8) :: ratio

  ! -------------------------------------

  smpm = sm_p -1


  DO sm_k=1,smpm

     ratio = DBLE(sm_k)/DBLE(sm_p)

     isa = 0
     DO is=1,nsp
        DO ia=1,na(is)
           isa = isa + 1
           DO i=1,3

              state(sm_k)%d3(i,isa) = state(0)%d3(i,isa)*(1.d0-ratio) &
                   & + state(sm_p)%d3(i,isa)*ratio
           ENDDO
        ENDDO
     ENDDO

  ENDDO

  RETURN

END SUBROUTINE LINEARP


!===============================================================
!
! subroutine ARC
!
!
!       This calculates the arclengths
!
!       key = 0 : arc from rep(0) 
!       key = 1 : diff in arc
!
!===============================================================

SUBROUTINE ARC(state,alpha,t_alpha,key)


  use ions_base, ONLY: na, nsp, nat
  USE path_variables, only: &
        sm_p => smd_p, &
        ptr => smd_ptr


  IMPLICIT NONE

  integer :: i,j,is,ia,sm_k,smpm, isa
  integer, intent(in) :: key

  type(ptr) :: state(0:sm_p)

  real(8), intent(out) :: alpha(0:sm_p),t_alpha 
  real(8) :: tmp(3,nat), dalpha(0:sm_p) 

  ! -------------------------------------------------


  smpm = sm_p -1

  dalpha(0:sm_p) = 0.d0
  alpha(0:sm_p) = 0.d0


  ! calculates arclength between replicas ------


  DO sm_k=1,sm_p

     isa = 0
     DO is=1,nsp
        DO ia=1,na(is)
           isa = isa + 1
           DO i=1,3

              tmp(i,isa) = state(sm_k)%d3(i,isa)-state(sm_k-1)%d3(i,isa)

           ENDDO
        ENDDO
     ENDDO

     isa = 0
     DO is=1,nsp
        DO ia=1,na(is)
           isa = isa + 1
           DO i=1,3

              dalpha(sm_k) = dalpha(sm_k) + tmp(i,isa)**2.d0

           ENDDO
        ENDDO
     ENDDO

     dalpha(sm_k) = SQRT(dalpha(sm_k))

  ENDDO


  ! calculates the arc length from rep(0)-----


  t_alpha = 0.d0

  alpha(0:sm_p) = 0.d0


  DO sm_k=0,sm_p
     DO j=0,sm_k
        alpha(sm_k) = alpha(sm_k) + dalpha(j)
     ENDDO
  ENDDO



  ! calculates total. arclength --------------

  t_alpha = alpha(sm_p)


  ! Normalize -------------------------------

  DO sm_k=1,sm_p
     dalpha(sm_k) = dalpha(sm_k)/t_alpha
     alpha(sm_k) = alpha(sm_k) / t_alpha
  ENDDO

  ! Choice -------------------

  IF(key==1) alpha(0:sm_p) = dalpha(0:sm_p)  



  RETURN


END SUBROUTINE ARC


!============================================================
!
! Subroutine ABSVEC
!
!  which calculates the absolute value |vec| of a vector 
!
!============================================================

SUBROUTINE ABSVEC(vec,norm)


  use ions_base, ONLY: na, nsp, nat

  IMPLICIT NONE

  integer :: i,is,ia,isa
  real(8),intent(in) :: vec(3,nat)
  real(8),intent(out) :: norm 


  norm = 0.d0

  isa = 0
  DO is=1,nsp
     DO ia=1,na(is)
        isa = isa + 1
        DO i=1,3

           norm = norm + vec(i,isa)**2.d0

        ENDDO
     ENDDO
  ENDDO

  norm = SQRT(norm)

  RETURN

END SUBROUTINE ABSVEC


!============================================================
!
! Subroutine Startcoord
!
! This subroutine is used to make the 6 degress of freedom
! the same. 
!
!
!============================================================

SUBROUTINE IMPOSECD(state,als,bes,chs)


  use ions_base, ONLY: na, nsp, nat
  use io_global, ONLY: stdout 

  IMPLICIT NONE

  integer :: i,is,ia
  integer,intent(in) :: als,bes,chs
  real(8) :: state(3,nat)
  real(8) :: dotp,theta,mag1,mag2


  write(stdout,*) " STARTCOORD used : ",als,bes,chs


  ! Translate --------------

  call tran(state,-state(1,als),'x')
  call tran(state,-state(2,als),'y')
  call tran(state,-state(3,als),'z')


  ! Rotate to 2(X,0,0) ---------

  ! to rotate around Z axis set y = 0
  ! find theta which is between (1,0,0) and XY proj of 2(x,y,z)

  mag1 = DSQRT(state(1,bes)**2.d0+state(2,bes)**2.d0)
  theta = ACOS(state(1,bes)/mag1)

  if(state(1,bes) >= 0.d0) then
     if(state(2,bes) >= 0.d0) then
        theta = -theta
     else
        theta = theta
     endif
  else
     if(state(2,bes) >= 0.d0) then
        theta = -theta
     else
        theta = theta
     endif
  endif


  call rot(state,theta,'z')


  ! to rotate around y set z = 0
  ! find theta which is between (1,0,0) and XZ proj of 2(x,0,z)


  mag1 = DSQRT(state(1,bes)**2.d0+state(3,bes)**2.d0)
  theta = ACOS(state(1,bes)/mag1)

  if(state(1,bes) >= 0.d0) then
     if(state(3,bes) >= 0.d0) then
        theta = theta
     else
        theta = -theta
     endif
  else
     if(state(3,bes) >= 0.d0) then
        theta = -theta
     else
        theta = theta
     endif
  endif

  call rot(state,theta,'y')



  ! Rotate to 3(X,Y,0) ---------

  ! to rotate around x set z = 0
  ! find theta which is between (0,1,0) and YZ proj of 3(x,y,z)

  mag1 = DSQRT(state(2,chs)**2.d0+state(3,chs)**2.d0)
  theta = ACOS(state(2,chs)/mag1)

  if(state(2,chs) >= 0.d0) then
     if(state(3,chs) >= 0.d0) then
        theta = -theta
     else
        theta = theta
     endif
  else
     if(state(3,chs) >= 0.d0) then
        theta = -theta
     else
        theta = theta
     endif
  endif

  call rot(state,theta,'x')

  RETURN

END SUBROUTINE IMPOSECD

!=========================================================
!
! Subroutine for rotations & translations 
!
! counter-clock wise
!
!=========================================================

subroutine rot(state,phi,direction)

  use ions_base, ONLY: na, nsp, nat
  use io_global, ONLY: stdout


  IMPLICIT NONE

  integer :: i,j,k,isa
  integer :: ia,is
  integer :: a,b,c
  real(8) :: rotm(3,3,3)
  real(8) :: state(3,nat)
  real(8), allocatable :: cd(:,:)
  real(8), allocatable :: newcd(:,:)
  real(8) :: tmp
  real(8), intent(in) :: phi
  real(8), parameter :: zero = 1.d-10
  character,intent(in) :: direction


  write(stdout,*) "ROT : ", direction

  allocate(cd(3,nat))
  allocate(newcd(3,nat))


  j=0
  DO is=1,nsp
     DO ia=1,na(is)
        j = j+1
        DO i=1,3
           cd(i,j) = state(i,j) 
        ENDDO
     ENDDO
  ENDDO


  !--- Building rotational matrix --------------!

  ! X -------

  k = 1

  rotm(1,1,k) = 1.d0
  rotm(2,1,k) = 0.d0
  rotm(3,1,k) = 0.d0

  rotm(1,2,k) = 0.d0
  rotm(2,2,k) = COS(phi)
  rotm(3,2,k) = SIN(phi)

  rotm(1,3,k) = 0.d0
  rotm(2,3,k) = -SIN(phi)
  rotm(3,3,k) = COS(phi)

  ! Y -------

  k = 2

  rotm(1,1,k) = COS(phi)
  rotm(2,1,k) = 0.d0
  rotm(3,1,k) = -SIN(phi)

  rotm(1,2,k) = 0.d0
  rotm(2,2,k) = 1.d0
  rotm(3,2,k) = 0.d0

  rotm(1,3,k) = SIN(phi)
  rotm(2,3,k) = 0.d0
  rotm(3,3,k) = COS(phi)


  ! Z --------

  k = 3

  rotm(1,1,k) = COS(phi)
  rotm(2,1,k) = SIN(phi)
  rotm(3,1,k) = 0.d0

  rotm(1,2,k) = -SIN(phi)
  rotm(2,2,k) = COS(phi)
  rotm(3,2,k) = 0.d0

  rotm(1,3,k) = 0.d0
  rotm(2,3,k) = 0.d0
  rotm(3,3,k) = 1.d0

  !--------------------------------!

  IF(direction == 'x') k = 1
  IF(direction == 'y') k = 2
  IF(direction == 'z') k = 3


  DO i=1,nat
     DO j=1,3
        tmp = rotm(j,1,k)*cd(1,i) + rotm(j,2,k)*cd(2,i) + rotm(j,3,k)*cd(3,i)
        if(dabs(tmp) < zero) then
           tmp = 0.d0
        endif
        newcd(j,i) = tmp
     ENDDO
  ENDDO

  j=0
  DO is=1,nsp
     DO ia=1,na(is)
        j = j+1
        DO i=1,3
           state(i,j) = newcd(i,j)
        ENDDO
     ENDDO
  ENDDO

  deallocate(cd)
  deallocate(newcd)

  RETURN

END SUBROUTINE ROT

!=====================================================

SUBROUTINE TRAN(state,dist,direction)

  use ions_base, ONLY: na, nsp, nat
  use io_global, ONLY: stdout

  IMPLICIT NONE

  integer :: i,j,k,is,ia,isa
  real(8) :: state(3,nat)
  real(8),intent(in) :: dist 
  character,intent(in) :: direction

  write(stdout,*) "TRAN : ", direction

  if(direction=='x') k=1
  if(direction=='y') k=2
  if(direction=='z') k=3

  isa = 0
  DO is=1,nsp
     DO ia=1,na(is)
        isa = isa + 1
        state(k,isa) = state(k,isa) + dist 
     ENDDO
  ENDDO

  RETURN

END SUBROUTINE TRAN

!=================================================
!
!  sub. distatoms - to calculate the closest atom dist.
!
!==================================================

SUBROUTINE DISTATOMS(state,mindist,isa1,isa2)

  use ions_base, ONLY: na, nsp, nat
  use cell_base, only: ainv, a1, a2, a3

  real(8),intent(in) :: state(3,nat)
  real(8),intent(out) :: mindist
  integer,intent(out) :: isa1, isa2
  real(8) :: dist, in(3), out(3) 
  integer :: is, ia, iis, iia, isa, iisa

  mindist  = 1.d10

  isa = 0
  DO is=1,nsp
     DO ia=1,na(is) 
        isa = isa + 1
        iisa = 0
        DO iis=is+1,nsp
           DO iia=1,na(iis)
              iisa = iisa + 1
              dist = 0.d0
              in(:) = state(:,isa)-state(:,iisa)
              call pbc(in,a1,a2,a3,ainv,out)
              dist = out(1)**2.d0 + out(2)**2.d0 + out(3)**2.d0 
              dist = SQRT(dist)
              IF(dist < mindist) THEN
                 mindist = dist 
                 isa1 = isa
                 isa2 = iisa
              ENDIF
           ENDDO
        ENDDO
     ENDDO
  ENDDO

  RETURN

END SUBROUTINE DISTATOMS


!===============================================================
!
!
subroutine init_path(sm_p,kwnp,stcd,nsp,nat,alat,nbeg,key) 

  USE input_parameters, only: sp_pos, pos, rd_pos, &
       & stcd1 => smd_stcd1, stcd2 => smd_stcd2, stcd3 => smd_stcd3, &
       & atomic_positions 

  USE parameters, ONLY: nsx
  USE smd_rep
  USE path_variables, only: &
        ptr => smd_ptr
  USE io_global, ONLY: stdout    
  USE cp_main_variables, ONLY: taub    

  implicit none


  INTEGER,intent(in) :: sm_p, kwnp,nsp, nat, nbeg, key 
  INTEGER :: sm_k, isa, na( nsx )
  INTEGER :: ia, is
  INTEGER :: als,bls,cls
  INTEGER :: i,j,k,redismi

  LOGICAL, intent(in) :: stcd

  TYPE(ptr), allocatable :: p_tau0(:) 
  real (8), allocatable :: guess(:,:,:)
  real (8), allocatable :: guess2(:,:,:)
  real (8), allocatable :: arc_in(:)
  real (8), allocatable :: arc_out(:)
  real (8), allocatable :: darc(:)

  real (8), intent(in) :: alat
  real (8) :: trashd

  ! key = 1  : SMOPT
  ! key = 2  : LINI
  ! key = 3  : POLM
  ! key = 4  : read all 


  ALLOCATE(rep(0:sm_p))
  ALLOCATE(rep_el(1:sm_p-1))

  DO sm_k=0,sm_p
     ALLOCATE(rep(sm_k)%tau0(3,nat))
     ALLOCATE(rep(sm_k)%taup(3,nat))
     ALLOCATE(rep(sm_k)%taum(3,nat))
     ALLOCATE(rep(sm_k)%taus(3,nat))
     ALLOCATE(rep(sm_k)%tausp(3,nat))
     ALLOCATE(rep(sm_k)%tausm(3,nat))
     ALLOCATE(rep(sm_k)%fion(3,nat))
     ALLOCATE(rep(sm_k)%fionm(3,nat))
     ALLOCATE(rep(sm_k)%fionp(3,nat))
     ALLOCATE(rep(sm_k)%vels(3,nat))
     ALLOCATE(rep(sm_k)%velsm(3,nat))
     ALLOCATE(rep(sm_k)%velsp(3,nat))
     ALLOCATE(rep(sm_k)%tan(3,nat))
  ENDDO

  IF( ALLOCATED( taub ) ) DEALLOCATE( taub )
  ALLOCATE(taub(3,nat))

  ! compute na(is)

  na = 0
  DO is = 1, nsp
    na( is ) = COUNT( sp_pos(1:nat) == is )
  END DO


  SELECT CASE ( key )

  CASE(1)  !  ... SMOPT --------


     IF(sm_p /= 3) call errore('path_smopt', ' sm_p /=3 with smopt', sm_p) 

     rep(0)%tau0 = 0.d0
     rep(3)%tau0 = 0.d0

     DO sm_k = 1, 2

        isa = 0
        DO is=1,nsp
           DO ia=1,nat

              rd_pos(:,ia) = pos(( 3 * ia - 2):( 3 * ia ),sm_k)

              IF( sp_pos(ia) == is ) THEN 

                 isa    = isa + 1
                 if( isa > nat) call errore(' init_path ',' isa > nat', isa )

                 rep(sm_k)%tau0(:,isa) = rd_pos(:, ia) 

                 IF(ia==stcd1) THEN
                    als = isa 
                 ELSEIF(ia==stcd2) THEN
                    bls = isa
                 ELSEIF(ia==stcd3) THEN
                    cls = isa
                 ENDIF

              ENDIF

           ENDDO
        ENDDO
     ENDDO

     IF(stcd) THEN
        CALL imposecd(rep(1)%tau0,als,bls,cls)
        CALL imposecd(rep(2)%tau0,als,bls,cls)
     ENDIF


  CASE(2) ! ... LINI ---------

     ALLOCATE(p_tau0(0:sm_p))


     ! ... pointer used later

     DO sm_k=0,sm_p 
        p_tau0(sm_k)%d3 => rep(sm_k)%tau0
     ENDDO


     isa = 0
     DO is=1,nsp
        DO ia=1,nat

           ! ... the 1st replica

           rd_pos(:,ia) = pos(( 3 * ia - 2):( 3 * ia ),1)

           IF( sp_pos(ia) == is ) THEN
              isa = isa + 1
              if( isa > nat) call errore(' init_path ',' isa > nat', isa )
              rep(0)%tau0(:,isa) = rd_pos(:, ia)
           ENDIF

           ! ... the last replica 

           rd_pos(:,ia) = pos(( 3 * ia - 2):( 3 * ia ),2)

           IF( sp_pos(ia) == is ) THEN

              rep(sm_p)%tau0(:,isa) = rd_pos(:, ia)

              IF(ia==stcd1) THEN
                 als = isa
              ELSEIF(ia==stcd2) THEN
                 bls = isa
              ELSEIF(ia==stcd3) THEN
                 cls = isa
              ENDIF

           ENDIF

        ENDDO
     ENDDO


     IF(stcd) THEN
        CALL imposecd(rep(0)%tau0,als,bls,cls)
        CALL imposecd(rep(sm_p)%tau0,als,bls,cls)
     ENDIF


     ! ... Linear interpolation

     CALL linearp(p_tau0)


     WRITE( stdout, *) "  " 
     WRITE( stdout, *) " INIT_PATH : linear interpolation used " 
     WRITE( stdout, *) "  " 

     ! ... Clean up !!

     DO sm_k=0,sm_p
        NULLIFY(p_tau0(sm_k)%d3)
     ENDDO

     DEALLOCATE(p_tau0)



  CASE(3)  ! ... POLM ---------


     ALLOCATE(p_tau0(0:sm_p))
     ALLOCATE(guess(3,nat,0:kwnp-1))
     ALLOCATE(guess2(3,nat,0:sm_p))
     ALLOCATE(arc_in(0:kwnp-1))
     ALLOCATE(arc_out(0:sm_p))
     ALLOCATE(darc(0:kwnp-1))

     ! ... pointer used later

     DO sm_k=0,sm_p
        p_tau0(sm_k)%d3 => rep(sm_k)%tau0
     ENDDO


     DO j=1,kwnp
        isa = 0
        DO is=1,nsp
           DO ia=1,nat

              rd_pos(:,ia) = pos(( 3 * ia - 2):( 3 * ia ),j)
              IF( sp_pos(ia) == is) THEN
                 isa = isa + 1
                 if( isa > nat) call errore(' init_path ',' isa > nat', isa )
                 guess(:,isa,j-1) = rd_pos(:, ia)
                 IF(ia==stcd1) THEN
                    als = isa
                 ELSEIF(ia==stcd2) THEN
                    bls = isa
                 ELSEIF(ia==stcd3) THEN
                    cls = isa
                 ENDIF

              ENDIF

           ENDDO
        ENDDO
     ENDDO


     ! ... stcd ...


     IF(stcd) THEN
        DO j=0,kwnp-1
           CALL imposecd(guess(:,:,j),als,bls,cls)
        ENDDO
     ENDIF


     ! ... Polm + grid ...

     ! ... Ideal arc


     DO sm_k=0,sm_p
        arc_out(sm_k) = DBLE(sm_k) * 1.d0/DBLE(sm_p)
     ENDDO

     ! ... Arc of kwon points


     darc(0) = 0.d0

     DO j=1,kwnp-1
        darc(j) = 0.d0

        isa = 0
        DO is=1,nsp
           DO ia=1,na(is)
              isa = isa + 1
              DO i=1,3
                 darc(j) = darc(j) + (guess(i,isa,j) - guess(i,isa,j-1))**2.d0
              ENDDO
           ENDDO
        ENDDO
        darc(j) = DSQRT(darc(j))
     ENDDO

     arc_in = 0.d0

     DO j=0,kwnp-1
        DO k=0,j
           arc_in(j) = arc_in(j) + darc(k)
        ENDDO
     ENDDO

     DO j=1,kwnp-1
        arc_in(j) = arc_in(j) /arc_in(kwnp-1)
     ENDDO

     ! ... Polynomial interpolation,  Neville's algorithm ... 


     DO sm_k=1,sm_p-1
        isa = 0
        DO is=1,nsp
           DO ia=1,na(is)
              isa = isa + 1
              DO i=1,3

                 call INTPOL_POLINT(arc_in(0:kwnp-1),guess(i,isa,0:kwnp-1), &
                      & kwnp, arc_out(sm_k), guess2(i,isa,sm_k),trashd)

                 rep(sm_k)%tau0(i,isa) = guess2(i,isa,sm_k)

              ENDDO
           ENDDO
        ENDDO
     ENDDO


     ! ...  initial and final

     isa = 0
     DO is=1,nsp
        DO ia=1,na(is)
           isa = isa + 1
           DO i=1,3
              rep(0)%tau0(i,isa) = guess(i,isa,0)
              rep(sm_p)%tau0(i,isa) = guess(i,isa,kwnp-1)
           ENDDO
        ENDDO
     ENDDO


     call ARC(p_tau0,arc_out,trashd,1)


     WRITE( stdout, *) " "
     WRITE( stdout, *) " INIT_PATH : polynomial interpolation used "
     WRITE( stdout, *) " "

     DO sm_k=0,sm_p
        WRITE( stdout, *) arc_out(sm_k)
     ENDDO


     ! ... Grid distribution ... 

     redismi = 10  ! Just a parameter for iteration

     DO j=1,redismi
        call REDIS(p_tau0)
     ENDDO

     call ARC(p_tau0,arc_out,trashd,1)

     WRITE( stdout, *) " "
     WRITE( stdout, *) " INIT_PATH : Grid distribution used "
     WRITE( stdout, *) " "

     DO sm_k=0,sm_p
        WRITE( stdout, *) arc_out(sm_k)
     ENDDO


     ! ... Clean up !!!

     DO sm_k=0,sm_p
        NULLIFY(p_tau0(sm_k)%d3)
     ENDDO

     DEALLOCATE(p_tau0)
     DEALLOCATE(guess)
     DEALLOCATE(guess2)
     DEALLOCATE(arc_in)
     DEALLOCATE(arc_out)
     DEALLOCATE(darc)

  CASE(4)  ! ... READ ALL ------

     nbeg_selection : IF(nbeg < 0) THEN 

        DO sm_k=0,sm_p
           isa = 0
           DO is=1,nsp
              DO ia=1,nat

                 rd_pos(:,ia) = pos(( 3 * ia - 2):( 3 * ia ),sm_k+1)
                 IF( sp_pos(ia) == is) THEN
                    isa = isa + 1
                    if( isa > nat) call errore(' init_path ',' isa > nat', isa )
                    rep(sm_k)%tau0(i,isa) = rd_pos(i, ia)
                    IF(ia==stcd1) THEN
                       als = isa
                    ELSEIF(ia==stcd2) THEN
                       bls = isa
                    ELSEIF(ia==stcd3) THEN
                       cls = isa
                    ENDIF
                 ENDIF

              ENDDO
           ENDDO
        ENDDO


        ! ... stcd ...


        IF(stcd) THEN
           DO sm_k=0,sm_p
              CALL imposecd(rep(sm_k)%tau0,als,bls,cls)
           ENDDO
        ENDIF


     ELSE

        isa = 0
        DO is=1,nsp
           DO ia=1,nat
              rd_pos(:,ia) = pos(( 3 * ia - 2):( 3 * ia ),1)
              IF( sp_pos(ia) == is) THEN
                 isa = isa + 1
                 if( isa > nat) call errore(' init_path ',' isa > nat', isa )
                 rep(0)%tau0(:,isa) = rd_pos(:, ia)
              ENDIF
              rd_pos(:,ia) = pos(( 3 * ia - 2):( 3 * ia ),sm_p+1)
              IF( sp_pos(ia) == is) THEN
                 rep(sm_p)%tau0(:,isa) = rd_pos(:, ia)
                 IF(ia==stcd1) THEN
                   als = isa
                 ELSEIF(ia==stcd2) THEN
                   bls = isa
                 ELSEIF(ia==stcd3) THEN
                   cls = isa
                 ENDIF
              ENDIF

           ENDDO
        ENDDO


        IF(stcd) THEN
           CALL imposecd(rep(0)%tau0,als,bls,cls)
           CALL imposecd(rep(sm_p)%tau0,als,bls,cls)
        ENDIF


     ENDIF nbeg_selection

  END SELECT


  ! ====== CONVERSION ========== 


  SELECT CASE ( atomic_positions )
     !
     !  convert input atomic positions to internally used format:
     !  tau0 in atomic units
     !
  CASE ('alat')
     !
     !  input atomic positions are divided by a0
     !
     DO sm_k=0,sm_p
        rep(sm_k)%tau0 = rep(sm_k)%tau0 * alat
     ENDDO
     !
  CASE ('bohr')
     !
     !  input atomic positions are in a.u.: do nothing
     !
     continue
  CASE ('crystal')
     !
     !  input atomic positions are in crystal axis ("scaled"):
     !
     CALL errore(' iosys ',' atomic_positions='//trim(atomic_positions)// &
          ' not implemented for SMD ', 1 )
     ! 
     !
  CASE ('angstrom')
     !
     !  atomic positions in A
     !
     DO sm_k=0,sm_p
        rep(sm_k)%tau0 = rep(sm_k)%tau0 / 0.529177
     ENDDO
     !
  CASE DEFAULT
     CALL errore(' iosys ',' atomic_positions='//trim(atomic_positions)// &
          ' not implemented ', 1 )
  END SELECT


  return

end subroutine init_path



subroutine allocate_path
  return
end subroutine allocate_path

