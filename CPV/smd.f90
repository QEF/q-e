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
! sub. r_to_s
!
!
!
!-----------------------------------------------------------------------
subroutine sminit (ibrav,celldm, ecut, ecutw,tranp,amprp,ndr,nbeg,  &
     tfirst,twmass,thdiag,iforceh,delt)
  !-----------------------------------------------------------------------
  !
  !     initialize G-vectors and related quantities
  !     use ibrav=0 for generic cell vectors given by the matrix h(3,3)
  !
  use control_flags, only: iprint, thdyn
  use io_global, only: stdout
  use gvec
  use gvecw, only: ngw
  use ions_base, only: na, pmass, nsp
  use cell_base, only: ainv, a1, a2, a3
  use elct
  use constants, only: pi, fpi
  use cell_base, only: wmass, hold, h
  use gvecw, only: ggp, agg => ecutz, sgg => ecsig, e0gg => ecfix
  use betax, only: mmx, refg
  !use restartsm
  use restart
  use parameters, only: nacx, nsx, natx
  USE smd_rep, only: rep
  USE smd_variables, only: sm_p

  implicit none
  ! input/output
  integer ibrav, ndr, nbeg, iforceh(3,3)
  logical tranp(nsx), tfirst, twmass, thdiag
  real(kind=8) amprp(nsx)
  real(kind=8) celldm(6), ecut, ecutw
  real(kind=8) delt
  ! local
  real(kind=8) randy
  integer i, j, ia, is, nfi


  ! YK
  ! present in the call to read(p)file, not actually used
  !      complex(kind=8) c0(1,1),cm(1,1)
  !      real(kind=8) taum(1,1,1),vel(1,1,1),velm(1,1,1),acc(nacx)
  !      real(kind=8) lambda(1,1),lambdam(1,1)
  !      real(kind=8) xnhe0,xnhem,vnhe,xnhp0,xnhpm,vnhp, ekincm
  !      real(kind=8) xnhh0(3,3),xnhhm(3,3),vnhh(3,3),velh(3,3)
  !      real(kind=8) fion(1,1,1)

  ! YK
  complex(kind=8) c0(ngw,nx),cm(ngw,nx)
  real(kind=8) taum(3,natx,nsx),vel(3,natx,nsx),velm(3,natx,nsx),acc(nacx)
  real(kind=8) lambda(nx,nx),lambdam(nx,nx)
  real(kind=8) xnhe0,xnhem,vnhe,xnhp0,xnhpm,vnhp, ekincm
  real(kind=8) xnhh0(3,3),xnhhm(3,3),vnhh(3,3),velh(3,3)
  real(kind=8) fion(3,natx,nsx)
  !


  integer :: sm_k,smpm 
  !
  !
  !     ==============================================================
  !     ==== generate reference g-space                           ==== 
  !     ==============================================================
  !
  smpm = sm_p -1
  !
  !
  !     ... initialize G-vectors and related quantities
  !
  call init1 ( rep(0)%tau0, ibrav, celldm, ecutw, ecut )

  !
  !
  ! taus = scaled, tau0 = alat units
  !

  DO sm_k=0, sm_p
     do is=1,nsp
        do ia=1,na(is)
           do i=1,3
              rep(sm_k)%taus(i,ia,is)=ainv(i,1)*rep(sm_k)%tau0(1,ia,is)       &
                   +ainv(i,2)*rep(sm_k)%tau0(2,ia,is)                 &
                   +ainv(i,3)*rep(sm_k)%tau0(3,ia,is)
           end do
        end do
     end do
  ENDDO
  !
  !
  refg=1.0*ecut/(mmx-1)
  WRITE( stdout,*) '   NOTA BENE: refg, mmx = ',refg,mmx
  !
  if(thdyn) then
     if(thdiag) then
        iforceh=0
        do i=1,3
           iforceh(i,i)=1
        enddo
     else
        iforceh=1
     endif
  endif
  !
  if( nbeg >= 0 ) then
     !
     ! read only h and hold from file ndr
     !
     call readfile_new                                              &
          &     (-1,ndr+1,h,hold,nfi,c0,cm,rep(0)%tau0,taum,vel,velm,acc,   &
          &       lambda,lambdam,xnhe0,xnhem,vnhe,xnhp0,xnhpm,vnhp,ekincm,   &
          &       xnhh0,xnhhm,vnhh,velh,ecut,ecutw,delt,pmass,ibrav,celldm,fion)
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
  allocate( ggp(ngw) )
  !
  !     ==============================================================
  !     ==== generate true g-space                                ==== 
  !     ==============================================================
  !
  call newinit( ibrav )
  !
  !
  DO sm_k = 0,sm_p
     !
     do is=1,nsp
        if(tranp(is)) then
           do ia=1,na(is)
              do i=1,3
                 rep(sm_k)%taus(i,ia,is)=rep(sm_k)%taus(i,ia,is)+amprp(is)*(randy()-0.5)
              end do
           end do
           !
           !     true tau (tau0) from scaled tau (taus)
           !
           do ia=1,na(is)
              do i=1,3
                 rep(sm_k)%tau0(i,ia,is) = h(i,1)*rep(sm_k)%taus(1,ia,is) &
                      + h(i,2)*rep(sm_k)%taus(2,ia,is) &
                      + h(i,3)*rep(sm_k)%taus(3,ia,is)
              end do
           end do
        end if
     end do
     !
  ENDDO
  !
  !
  if(.not. twmass) then
     WRITE( stdout,998) wmass
  else
     wmass=0.
     do is=1,nsp
        wmass=wmass+na(is)*pmass(is)
     enddo
     wmass=wmass*0.75/pi/pi
     WRITE( stdout,999) wmass
  endif
998 format(' wmass (read from input) = ',f15.2,/)
999 format(' wmass (calculated) = ',f15.2,/)
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


  use ions_base, ONLY: na, nsp
  use smd_ene, ONLY: etot_ar
  use smd_variables

  IMPLICIT NONE

  integer :: i,j,is,ia,sm_k,smpm

  type(ptr) :: state(0:sm_p)
  type(ptr) :: tan(0:sm_p)

  real(kind=8) :: ene1,ene2,tmp1,tmp2
  real(kind=8), parameter :: zero = 1.d-15

  ! ---------------------------- !


  tan(0)%d3(:,:,:) = 0.d0
  tan(sm_p)%d3(:,:,:) = 0.d0   

  smpm = sm_p -1


  DO sm_k=1,smpm      ! >>>>>>>>>>>>>>>>>>>>>>  !


     IF((etot_ar(sm_k+1) >= etot_ar(sm_k)) .AND. &
          &  (etot_ar(sm_k) >= etot_ar(sm_k-1))) THEN

        DO is=1,nsp
           DO ia=1,na(is)
              DO i=1,3
                 tan(sm_k)%d3(i,ia,is)=state(sm_k+1)%d3(i,ia,is) & 
                      & -state(sm_k)%d3(i,ia,is)
              ENDDO
           ENDDO
        ENDDO

     ELSE IF((etot_ar(sm_k+1) <= etot_ar(sm_k)) .AND. &
          & (etot_ar(sm_k) <= etot_ar(sm_k-1))) THEN

        DO is=1,nsp
           DO ia=1,na(is)
              DO i=1,3
                 tan(sm_k)%d3(i,ia,is)=state(sm_k)%d3(i,ia,is) &
                      & -state(sm_k-1)%d3(i,ia,is)
              ENDDO
           ENDDO
        ENDDO

     ELSE

        tmp1=ABS(etot_ar(sm_k-1)-etot_ar(sm_k))
        tmp2=ABS(etot_ar(sm_k)-etot_ar(sm_k+1))
        ene1=MAX(tmp1,tmp2)
        ene2=MIN(tmp1,tmp2)

        IF(etot_ar(sm_k+1) >= etot_ar(sm_k-1)) THEN

           DO is=1,nsp
              DO ia=1,na(is)
                 DO i=1,3

                    tan(sm_k)%d3(i,ia,is) = ene1*(state(sm_k+1)%d3(i,ia,is) & 
                         & -state(sm_k)%d3(i,ia,is)) &
                         & + ene2*(state(sm_k)%d3(i,ia,is)  &
                         & -state(sm_k-1)%d3(i,ia,is))

                 ENDDO
              ENDDO
           ENDDO

        ELSE

           DO is=1,nsp
              DO ia=1,na(is)
                 DO i=1,3

                    tan(sm_k)%d3(i,ia,is) = ene2*(state(sm_k+1)%d3(i,ia,is) &
                         & -state(sm_k)%d3(i,ia,is)) &
                         & + ene1*(state(sm_k)%d3(i,ia,is)  &
                         & -state(sm_k-1)%d3(i,ia,is))

                 ENDDO
              ENDDO
           ENDDO

        ENDIF

     ENDIF


     ! Normalization of tangent ------------------------------- 

     tmp1=0.d0

     DO is=1,nsp
        DO ia=1,na(is)
           DO i=1,3

              tmp1 = tmp1 + tan(sm_k)%d3(i,ia,is)**2.d0

           ENDDO
        ENDDO
     ENDDO


     tmp1 = SQRT(tmp1)

     DO is=1,nsp
        DO ia=1,na(is)
           DO i=1,3

              IF(ABS(tan(sm_k)%d3(i,ia,is)) >= zero ) THEN
                 tan(sm_k)%d3(i,ia,is) = tan(sm_k)%d3(i,ia,is)/tmp1
              ELSE
                 tan(sm_k)%d3(i,ia,is)=0.d0
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


  use ions_base, ONLY: na, nsp
  use parameters, only: nsx,natx
  use smd_variables

  IMPLICIT NONE

  integer :: i,is,ia

  real(kind=8) :: vec(3,natx,nsx),tan(3,natx,nsx)
  real(kind=8) :: paraforce
  real(kind=8) :: dotp

  !----------------------------------------------


  ! calculate DOT product of the vectors -----

  dotp = 0.d0

  DO is=1,nsp
     DO ia=1,na(is)
        DO i=1,3
           dotp = dotp + vec(i,ia,is)*tan(i,ia,is)
        ENDDO
     ENDDO
  ENDDO

  paraforce = dotp


  ! calculate PERP of the vec  -----------

  DO is=1,nsp
     DO ia=1,na(is)
        DO i=1,3
           vec(i,ia,is) = vec(i,ia,is) - dotp*tan(i,ia,is)
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
  use parameters, only: nsx,natx
  use smd_variables


  IMPLICIT NONE

  integer :: i,j,is,ia,sm_k,smpm

  type(ptr) :: state(0:sm_p)

  real(kind=8) :: ratio

  ! -------------------------------------

  smpm = sm_p -1


  DO sm_k=1,smpm

     ratio = DBLE(sm_k)/DBLE(sm_p)

     DO is=1,nsp
        DO ia=1,na(is)
           DO i=1,3

              state(sm_k)%d3(i,ia,is) = state(0)%d3(i,ia,is)*(1.d0-ratio) &
                   & + state(sm_p)%d3(i,ia,is)*ratio
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


  use ions_base, ONLY: na, nsp
  use parameters, only: nsx,natx
  use smd_variables


  IMPLICIT NONE

  integer :: i,j,is,ia,sm_k,smpm
  integer, intent(in) :: key

  type(ptr) :: state(0:sm_p)

  real(kind=8), intent(out) :: alpha(0:sm_p),t_alpha 
  real(kind=8) :: tmp(3,natx,nsx), dalpha(0:sm_p) 

  ! -------------------------------------------------


  smpm = sm_p -1

  dalpha(0:sm_p) = 0.d0
  alpha(0:sm_p) = 0.d0


  ! calculates arclength between replicas ------


  DO sm_k=1,sm_p

     DO is=1,nsp
        DO ia=1,na(is)
           DO i=1,3

              tmp(i,ia,is) = state(sm_k)%d3(i,ia,is)-state(sm_k-1)%d3(i,ia,is)

           ENDDO
        ENDDO
     ENDDO

     DO is=1,nsp
        DO ia=1,na(is)
           DO i=1,3

              dalpha(sm_k) = dalpha(sm_k) + tmp(i,ia,is)**2.d0

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


  use ions_base, ONLY: na, nsp
  use parameters, only: nsx,natx
  use smd_variables

  IMPLICIT NONE

  integer :: i,is,ia
  real(kind=8),intent(in) :: vec(3,natx,nsx)
  real(kind=8),intent(out) :: norm 


  norm = 0.d0

  DO is=1,nsp
     DO ia=1,na(is)
        DO i=1,3

           norm = norm + vec(i,ia,is)**2.d0

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

SUBROUTINE IMPOSECD(state,al,als,be,bes,ch,chs)


  use ions_base, ONLY: na, nsp
  use parameters, only: nsx,natx
  use io_global, ONLY: stdout 

  IMPLICIT NONE

  integer :: i,is,ia
  integer,intent(in) :: al,be,ch,als,bes,chs
  real(kind=8) :: state(3,natx,nsx)
  real(kind=8) :: dotp,theta,mag1,mag2


  write(stdout,*) " STARTCOORD used : ",al,als,be,bes,ch,chs


  ! Translate --------------

  call tran(state,-state(1,al,als),'x')
  call tran(state,-state(2,al,als),'y')
  call tran(state,-state(3,al,als),'z')


  ! Rotate to 2(X,0,0) ---------

  ! to rotate around Z axis set y = 0
  ! find theta which is between (1,0,0) and XY proj of 2(x,y,z)

  mag1 = DSQRT(state(1,be,bes)**2.d0+state(2,be,bes)**2.d0)
  theta = ACOS(state(1,be,bes)/mag1)

  if(state(1,be,bes) >= 0.d0) then
     if(state(2,be,bes) >= 0.d0) then
        theta = -theta
     else
        theta = theta
     endif
  else
     if(state(2,be,bes) >= 0.d0) then
        theta = -theta
     else
        theta = theta
     endif
  endif


  call rot(state,theta,'z')


  ! to rotate around y set z = 0
  ! find theta which is between (1,0,0) and XZ proj of 2(x,0,z)


  mag1 = DSQRT(state(1,be,bes)**2.d0+state(3,be,bes)**2.d0)
  theta = ACOS(state(1,be,bes)/mag1)

  if(state(1,be,bes) >= 0.d0) then
     if(state(3,be,bes) >= 0.d0) then
        theta = theta
     else
        theta = -theta
     endif
  else
     if(state(3,be,bes) >= 0.d0) then
        theta = -theta
     else
        theta = theta
     endif
  endif

  call rot(state,theta,'y')



  ! Rotate to 3(X,Y,0) ---------

  ! to rotate around x set z = 0
  ! find theta which is between (0,1,0) and YZ proj of 3(x,y,z)

  mag1 = DSQRT(state(2,ch,chs)**2.d0+state(3,ch,chs)**2.d0)
  theta = ACOS(state(2,ch,chs)/mag1)

  if(state(2,ch,chs) >= 0.d0) then
     if(state(3,ch,chs) >= 0.d0) then
        theta = -theta
     else
        theta = theta
     endif
  else
     if(state(3,ch,chs) >= 0.d0) then
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

  use ions_base, ONLY: na, nsp
  use parameters, only: nsx,natx
  use io_global, ONLY: stdout


  IMPLICIT NONE

  integer :: i,j,k
  integer :: ia,is,nat
  integer :: a,b,c
  real(kind=8) :: rotm(3,3,3)
  real(kind=8) :: state(3,natx,nsx)
  real(kind=8), allocatable :: cd(:,:)
  real(kind=8), allocatable :: newcd(:,:)
  real(kind=8) :: tmp
  real(kind=8), intent(in) :: phi
  real(kind=8), parameter :: zero = 1.d-10
  character,intent(in) :: direction


  write(stdout,*) "ROT : ", direction

  nat=0
  DO is=1,nsp
     nat=nat+na(is)
  ENDDO

  allocate(cd(3,nat))
  allocate(newcd(3,nat))


  j=0
  DO is=1,nsp
     DO ia=1,na(is)
        j = j+1
        DO i=1,3
           cd(i,j) = state(i,ia,is) 
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
           state(i,ia,is) = newcd(i,j)
        ENDDO
     ENDDO
  ENDDO

  deallocate(cd)
  deallocate(newcd)

  RETURN

END SUBROUTINE ROT

!=====================================================

SUBROUTINE TRAN(state,dist,direction)

  use ions_base, ONLY: na, nsp
  use parameters, only: nsx,natx
  use io_global, ONLY: stdout

  IMPLICIT NONE

  integer :: i,j,k,is,ia
  real(kind=8) :: state(3,natx,nsx)
  real(kind=8),intent(in) :: dist 
  character,intent(in) :: direction

  write(stdout,*) "TRAN : ", direction

  if(direction=='x') k=1
  if(direction=='y') k=2
  if(direction=='z') k=3

  DO is=1,nsp
     DO ia=1,na(is)
        state(k,ia,is) = state(k,ia,is) + dist 
     ENDDO
  ENDDO

  RETURN

END SUBROUTINE TRAN

!=================================================
!
!  sub. distatoms - to calculate the closest atom dist.
!
!==================================================

SUBROUTINE DISTATOMS(state,mindist,is1,ia1,is2,ia2)

  use ions_base, ONLY: na, nsp
  use parameters, only: nsx,natx
  use cell_base, only: ainv, a1, a2, a3

  integer :: i,is,ia
  real(kind=8),intent(in) :: state(3,natx,nsx)
  real(kind=8),intent(out) :: mindist
  real(kind=8) :: dist, in(3), out(3) 
  integer,intent(out) :: is1,ia1,is2,ia2

  mindist  = 1.d10

  DO is=1,nsp
     DO ia=1,na(is) 
        DO iis=is+1,nsp
           DO iia=1,na(iis)
              dist = 0.d0
              DO i=1,3
                 in(i) = state(i,ia,is)-state(i,iia,iis)
              ENDDO
              call pbc(in,a1,a2,a3,ainv,out)
              dist = out(1)**2.d0 + out(2)**2.d0 + out(3)**2.d0 
              dist = SQRT(dist)
              IF(dist < mindist) THEN
                 mindist = dist 
                 is1 = is
                 ia1 = ia
                 is2 = iis
                 ia2 = iia
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

  USE parameters, ONLY: nsx, natx
  USE smd_rep
  USE smd_variables, ONLY: ptr
  USE io_global, ONLY: stdout    

  implicit none


  INTEGER,intent(in) :: sm_p, kwnp,nsp, nat, nbeg, key 
  INTEGER :: sm_k, isa, na( nsx )
  INTEGER :: index(2,3), ai, ia, is
  INTEGER :: al,als,bl,bls,cl,cls
  INTEGER :: i,j,k,redismi

  LOGICAL, intent(in) :: stcd

  TYPE(ptr), allocatable :: p_tau0(:) 
  real (kind=8), allocatable :: guess(:,:,:,:)
  real (kind=8), allocatable :: guess2(:,:,:,:)
  real (kind=8), allocatable :: arc_in(:)
  real (kind=8), allocatable :: arc_out(:)
  real (kind=8), allocatable :: darc(:)

  real (kind=8), intent(in) :: alat
  real (kind=8) :: trashd

  ! key = 1  : SMOPT
  ! key = 2  : LINI
  ! key = 3  : POLM
  ! key = 4  : read all 



  ALLOCATE(rep(0:sm_p))
  ALLOCATE(rep_el(1:sm_p-1))



  SELECT CASE ( key )



  CASE(1)  !  ... SMOPT --------


     IF(sm_p /= 3) call errore('path_smopt', ' sm_p /=3 with smopt', sm_p) 

     na = 0
     isa = 0
     ai = 0

     rep(0)%tau0 = 0.d0
     rep(3)%tau0 = 0.d0

     DO sm_k=1,2
        DO is=1,nsp
           DO ia=1,nat
              ai = ai +1

              IF(ai==stcd1) THEN
                 al = ia
                 als = is 
              ELSEIF(ai==stcd2) THEN
                 bl = ia
                 bls = is
              ELSEIF(ai==stcd3) THEN
                 cl = ia
                 cls = is
              ENDIF

              rd_pos(:,ia) = pos(( 3 * ia - 2):( 3 * ia ),sm_k)
              IF( sp_pos(ia) == is) THEN 
                 na(is) = na(is) + 1 
                 if( na(is) > natx) call errore(' cards',' na > natx', na(is) )
                 DO i=1,3
                    rep(sm_k)%tau0(i,na(is),is) = rd_pos(i, ia) 
                 ENDDO
              ENDIF

           ENDDO
        ENDDO
     ENDDO

     IF(stcd) THEN
        CALL imposecd(rep(1)%tau0,al,als,bl,bls,cl,cls)
        CALL imposecd(rep(2)%tau0,al,als,bl,bls,cl,cls)
     ENDIF


  CASE(2) ! ... LINI ---------

     ALLOCATE(p_tau0(0:sm_p))


     ! ... pointer used later

     DO sm_k=0,sm_p 
        p_tau0(sm_k)%d3 => rep(sm_k)%tau0
     ENDDO


     na = 0
     isa = 0
     ai = 0


     DO is=1,nsp
        DO ia=1,nat
           ai = ai +1

           IF(ai==stcd1) THEN
              al = ia
              als = is
           ELSEIF(ai==stcd2) THEN
              bl = ia
              bls = is
           ELSEIF(ai==stcd3) THEN
              cl = ia
              cls = is
           ENDIF


           ! ... the 1st replica

           rd_pos(:,ia) = pos(( 3 * ia - 2):( 3 * ia ),1)
           IF( sp_pos(ia) == is) THEN
              na(is) = na(is) + 1
              if( na(is) > natx) call errore(' cards',' na > natx', na(is) )
              DO i=1,3
                 rep(0)%tau0(i,na(is),is) = rd_pos(i, ia)
              ENDDO
           ENDIF


           ! ... the last replica 

           rd_pos(:,ia) = pos(( 3 * ia - 2):( 3 * ia ),2)
           IF( sp_pos(ia) == is) THEN
              if( na(is) > natx) call errore(' cards',' na > natx', na(is) )
              DO i=1,3
                 rep(sm_p)%tau0(i,na(is),is) = rd_pos(i, ia)
              ENDDO
           ENDIF

        ENDDO
     ENDDO


     IF(stcd) THEN
        CALL imposecd(rep(0)%tau0,al,als,bl,bls,cl,cls)
        CALL imposecd(rep(sm_p)%tau0,al,als,bl,bls,cl,cls)
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
     ALLOCATE(guess(3,natx,nsx,0:kwnp-1))
     ALLOCATE(guess2(3,natx,nsx,0:sm_p))
     ALLOCATE(arc_in(0:kwnp-1))
     ALLOCATE(arc_out(0:sm_p))
     ALLOCATE(darc(0:kwnp-1))

     ! ... pointer used later

     DO sm_k=0,sm_p
        p_tau0(sm_k)%d3 => rep(sm_k)%tau0
     ENDDO



     na = 0
     isa = 0

     DO j=1,kwnp
        ai = 0
        DO is=1,nsp
           DO ia=1,nat
              ai = ai +1

              IF(ai==stcd1) THEN
                 al = ia
                 als = is
              ELSEIF(ai==stcd2) THEN
                 bl = ia
                 bls = is
              ELSEIF(ai==stcd3) THEN
                 cl = ia
                 cls = is
              ENDIF

              rd_pos(:,ia) = pos(( 3 * ia - 2):( 3 * ia ),j)
              IF( sp_pos(ia) == is) THEN
                 na(is) = na(is) + 1
                 if( na(is) > natx) call errore(' cards',' na > natx', na(is) )
                 DO i=1,3
                    guess(i,na(is),is,j-1) = rd_pos(i, ia)
                 ENDDO
              ENDIF

           ENDDO
        ENDDO
     ENDDO


     ! ... stcd ...


     IF(stcd) THEN
        DO j=0,kwnp-1
           CALL imposecd(guess(:,:,:,j),al,als,bl,bls,cl,cls)
        ENDDO
     ENDIF


     ! ... Polm + grid ...

     ! ... Ideal arc


     DO sm_k=0,sm_p
        arc_out(sm_k) = dble(sm_k) * 1.d0/dble(sm_p)
     ENDDO

     ! ... Arc of kwon points


     darc(0) = 0.d0

     DO j=1,kwnp-1
        darc(j) = 0.d0

        DO is=1,nsp
           DO ia=1,na(is)
              DO i=1,3
                 darc(j) = darc(j) + (guess(i,ia,is,j) - guess(i,ia,is,j-1))**2.d0
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
        DO is=1,nsp
           DO ia=1,na(is)
              DO i=1,3

                 call INTPOL_POLINT(arc_in(0:kwnp-1),guess(i,ia,is,0:kwnp-1), &
                      & kwnp, arc_out(sm_k), guess2(i,ia,is,sm_k),trashd)

                 rep(sm_k)%tau0(i,ia,is) = guess2(i,ia,is,sm_k)

              ENDDO
           ENDDO
        ENDDO
     ENDDO


     ! ...  initial and final

     DO is=1,nsp
        DO ia=1,na(is)
           DO i=1,3
              rep(0)%tau0(i,ia,is) = guess(i,ia,is,0)
              rep(sm_p)%tau0(i,ia,is) = guess(i,ia,is,kwnp-1)
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
           DO is=1,nsp
              DO ia=1,nat
                 ai = ai +1

                 IF(ai==stcd1) THEN
                    al = ia
                    als = is
                 ELSEIF(ai==stcd2) THEN
                    bl = ia
                    bls = is
                 ELSEIF(ai==stcd3) THEN
                    cl = ia
                    cls = is
                 ENDIF

                 rd_pos(:,ia) = pos(( 3 * ia - 2):( 3 * ia ),sm_k+1)
                 IF( sp_pos(ia) == is) THEN
                    na(is) = na(is) + 1
                    if( na(is) > natx) call errore(' cards',' na > natx', na(is) )
                    DO i=1,3
                       rep(sm_k)%tau0(i,na(is),is) = rd_pos(i, ia)
                    ENDDO
                 ENDIF

              ENDDO
           ENDDO
        ENDDO


        ! ... stcd ...


        IF(stcd) THEN
           DO sm_k=0,sm_p
              CALL imposecd(rep(sm_k)%tau0,al,als,bl,bls,cl,cls)
           ENDDO
        ENDIF


     ELSE

        DO is=1,nsp
           DO ia=1,nat
              ai = ai +1

              IF(ai==stcd1) THEN
                 al = ia
                 als = is
              ELSEIF(ai==stcd2) THEN
                 bl = ia
                 bls = is
              ELSEIF(ai==stcd3) THEN
                 cl = ia
                 cls = is
              ENDIF

              rd_pos(:,ia) = pos(( 3 * ia - 2):( 3 * ia ),1)
              IF( sp_pos(ia) == is) THEN
                 na(is) = na(is) + 1
                 if( na(is) > natx) call errore(' cards',' na > natx', na(is) )
                 DO i=1,3
                    rep(0)%tau0(i,na(is),is) = rd_pos(i, ia)
                 ENDDO
              ENDIF

              rd_pos(:,ia) = pos(( 3 * ia - 2):( 3 * ia ),sm_p+1)
              IF( sp_pos(ia) == is) THEN
                 if( na(is) > natx) call errore(' cards',' na > natx', na(is) )
                 DO i=1,3
                    rep(sm_p)%tau0(i,na(is),is) = rd_pos(i, ia)
                 ENDDO
              ENDIF

           ENDDO
        ENDDO


        IF(stcd) THEN
           CALL imposecd(rep(0)%tau0,al,als,bl,bls,cl,cls)
           CALL imposecd(rep(sm_p)%tau0,al,als,bl,bls,cl,cls)
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


!===============================================================================!
! 
!

subroutine r_to_s(tau,taus)

  use cell_base, only: ainv
  use ions_base, only: na, nsp
  use parameters, only: natx, nsx

  implicit none

  integer :: is,ia,i
  real (kind=8), intent(out) :: taus(3,natx,nsx)
  real (kind=8), intent(in)  :: tau(3,natx,nsx)

  do is=1,nsp
     do ia=1,na(is)
        do i=1,3
           taus(i,ia,is)=ainv(i,1)*tau(1,ia,is)       &
                +ainv(i,2)*tau(2,ia,is)       &
                +ainv(i,3)*tau(3,ia,is)
        end do
     end do
  end do

  return

end subroutine r_to_s


!===============================================================================
!
!
