MODULE zero_mod


  USE kinds, ONLY: DP
  use splines
  use fft_base
  use becmod


  SAVE
!variabili per potenziale non locale
  logical  ::l_zero,l_non_loc,l_tab
  real(DP), allocatable :: tabr(:,:,:,:)
  real(DP), allocatable :: tabr_d2y(:,:,:,:)
  real(DP), allocatable :: tablocal_hg(:,:,:)
  real(DP), allocatable :: tablocal_d2y_hg(:,:,:)
  TYPE (bec_type) :: becpr(3)
  TYPE (bec_type) :: becpd(3)
  TYPE (bec_type) :: becprd(3,3)
  complex(DP), allocatable :: xvkb(:,:,:)
  complex(DP), allocatable :: dvkb(:,:,:)
  complex(DP), allocatable :: xdvkb(:,:,:,:)


  !the component of the current here computed
  real(dp) ::z_current(3), i_current(3), i_current_a(3), i_current_b(3), i_current_c(3), i_current_d(3), i_current_e(3)
 
  !variables depending on the step
  !wavefunction
  complex(DP),allocatable :: evc_uno(:,:)
  real(DP), allocatable   :: charge(:)

  !ion positions and velocities
  real(DP), allocatable ::ion_pos(:,:) 
  real(DP), allocatable ::ion_vel(:,:) 
  !second ion positions and velocities read from input
  real(DP), allocatable ::ion_pos2(:,:) ! must call convert_tau from ../PW/src/input.f90 to obtain correct units for positions
  real(DP), allocatable ::ion_vel2(:,:)
  character(len=256) :: second_vel_pos_fname

  !input from stdout
  integer        :: natoms !cutoff per somme in griglia reale
  integer        :: n_max !cutoff per somme in griglia reale
  character(256) :: status !what to do with the program. "initialize" or "compute"
  real(kind=DP)  :: eta !ewald factor for convergence

  
  real(DP) ::I_primo_rec
  real(DP), allocatable   :: H_g(:,:,:,:)          
  complex(DP), allocatable   :: u_g(:,:) !ngm,a    
  complex(DP), allocatable ::charge_g(:)
  real(DP), allocatable   :: I_uno_g(:,:,:) !griglia,a,b     
  real(DP), allocatable   :: I_due_g(:)  
  real(DP)  ::I_primo
 
  !variable for S_due
  real(DP), allocatable   ::alpha_0_all(:,:),alpha_2_all(:,:)


  !counters and indexes
  integer      ::nr3s_start, nr3s_end !z-planes in processor
  integer      ::nr3s_start_par
  contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real(kind=DP) function h_(x)
  real(kind=DP) :: x
  real(kind=DP), external :: qe_erfc
  h_=qe_erfc(x)/x
  end function h_

  real(kind=DP) function hp_(x)
  use constants, only :pi
  real(kind=DP) :: x
  real(kind=DP), external :: qe_erfc
  hp_=-(2.d0/sqrt(pi))*(1.d0/x)*exp(-x*x)-1.d0/(x*x)*qe_erfc(x)
  end function hp_


!!!!!!!!!!!!!!!!!!!!!!!!!!  

  real(kind=DP) function modulus(vector)
  real(kind=DP) ::vector(3)
  modulus=sqrt(vector(1)**2+vector(2)**2+vector(3)**2)
  end function modulus


  
!!!!!!!!!!!!!!!!!!!!!!!!



  real(DP) function alpha_0_lr(eta,G_sq,flag)
  real(DP)  ::eta,G_sq
  integer ::flag
  if (flag==1) then
      alpha_0_lr=exp(-G_sq/(4.d0*eta))/G_sq
  else
      alpha_0_lr=1.d0/G_sq
  end if
  end function alpha_0_lr



!!!!!!!!!!!!!!!!!!!!!!!



  real(DP) function alpha_2_lr(eta,G_sq,flag)
  use constants, only :pi
  real(kind=DP), external :: qe_erf
  real(DP)  ::eta,G_sq
  integer ::flag
  if (flag==1) then
      alpha_2_lr=1.d0/G_sq * ( (3.d0*sqrt(pi)*sqrt(eta/G_sq)) * qe_erf(sqrt(G_sq/(eta))/2.d0) - exp(-(G_sq)/(4.d0*eta)) )
  else
      alpha_2_lr=2.d0/G_sq
  end if
  end function alpha_2_lr



!!!!!!!!!!!!


  subroutine mid(vin,vout)
  use kinds,     only:DP
  use cell_base, only:alat
! apply the minimum image distance in a cubic cell
  implicit none
  real(DP), intent(in) :: vin(3)
  real(DP), intent(out) :: vout(3)
  !local
  integer :: i
  do i=1,3
    vout(i)=vin(i)-nint(vin(i)/alat)*alat
  end do
  end subroutine mid



!!!!!!!!!!!!!!!!!!!!!!!!!!!



  subroutine pbc(vin,vout)
  use kinds,     only:DP
  use cell_base, only:alat
! apply the minimum image distance in a cubic cell
  implicit none
  real(DP), intent(in) :: vin(3)
  real(DP), intent(out) :: vout(3)
  !local
  integer :: i,n
  vout(:)=0.d0
  do i=1,3
     if (vin(i)>=0) then
         n=int(vin(i)/alat)
     else
         n=int(vin(i)/alat)-1
     end if
     vout(i)=vin(i)-dble(n)*alat
  end do
  end subroutine pbc

subroutine check_positions(ion_pos)
  use ions_base, only :tau,nat
  use cell_base, only :alat
  real(DP), intent(in) ::ion_pos(3,nat)
  integer ::coord,iatom
  do iatom=1,nat
     do coord=1,3
        if (abs(tau(coord,iatom)*alat-ion_pos(coord,iatom))>1.E-4) then
             call errore('check_positions','positions from MD and from PW not matching',1)         
        end if
     end do
  end do   
end subroutine check_positions

!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!
  
  subroutine I_uno_value(y,x,a,b,flag)
  use cell_base, only:tpiba
  use gvect, only: ngm,gstart,g,gl,gg
  use mp,        only :mp_sum
  use wvfct, only :npw
  use mp_pools, only : intra_pool_comm
  implicit none
  integer, intent(in) ::flag
  real(DP), intent(out) ::y
  real(DP),intent(in)  ::x(3)
  integer,intent(in)   ::a,b
  integer ::igm

!  integer ::ng_max

  real(DP) :: comp_iso
  
  call start_clock( 'rec' )
!  ng_max=10

  y=0.d0
! do igm=gstart,ngm
  do igm=gstart,npw
     if	  ( (gg(igm)*tpiba*tpiba)/(4.d0*eta) > 20.d0)	exit
     y=y+2.d0*(I_uno_g(igm,a,b)*cos(DOT_PRODUCT(g(1:3,igm), x(1:3)) *tpiba) )
!     if (gl(igm)>ng_max) exit
  end do
  if (gstart==2) y=y+I_uno_g(1,a,b)
  call mp_sum(y, intra_pool_comm)
  call stop_clock( 'rec' )
  if (flag==1) then
      call add_local_uno(y,x,a,b)
  end if
!change for opt   
!  call I_due_value(comp_iso,x,1)
!  if (a==b) then
!     y=y-comp_iso
!  end if
   if (a==b) then
      call I_due_value(comp_iso,x,1)
      y=y-comp_iso
   end if 
 
  end subroutine I_uno_value

  subroutine add_local_uno(value,pos,a,b)
  use cell_base, only:at,alat
  use mp_world, only: nproc, mpime
  use mp, only:mp_sum
  use mp_pools, only : intra_pool_comm
  implicit none
  integer,  intent(in) ::a,b
  real(DP), intent(in) ::pos(3)
  real(DP), intent(inout) ::value
  real(DP) ::u(3),u_mod,n(3)
  integer  ::n_x,n_y,n_z
  integer   :: l_blk,nbegin,nend
  real(DP)  :: value1
!!!!!!!!!!!!!!!!!!!!!
  call start_clock( 'real' )
  l_blk= (2*n_max+1)/nproc
  if(l_blk*nproc <  (2*n_max+1) ) l_blk = l_blk+1
  nbegin=mpime*l_blk-n_max
  nend=nbegin+l_blk-1
  if(nend > n_max) nend = n_max
  value1=0.d0
  do n_x=nbegin,nend
     do n_y=-n_max,n_max
        do n_z=-n_max,n_max
           n(1:3) = n_x * at(1:3,1)*alat + n_y * at(1:3,2)*alat + n_z * at(1:3,3)*alat
           u(1:3) = pos(1:3)-n(1:3)
           u_mod=modulus(u)
           value1=value1+eta*hp_(sqrt(eta)*u_mod)*u(a)*u(b)/u_mod
           if (a==b) then
              value1=value1+sqrt(eta)*h_(sqrt(eta)*u_mod)
           end if
        end do
     end do
  end do
  call mp_sum(value1, intra_pool_comm)
  value=value+value1
  call stop_clock( 'real' )
!!!!!!!!!!!!!!!!!!! 
  !call start_clock( 'real' )
  !do n_x=-n_max,n_max  
  !   do n_y=-n_max,n_max 
  !      do n_z=-n_max,n_max 
  !         n(1:3) = n_x * at(1:3,1)*alat + n_y * at(1:3,2)*alat + n_z * at(1:3,3)*alat
  !         u(1:3) = pos(1:3)-n(1:3)
  !         u_mod=modulus(u) 
  !         value=value+eta*hp_(sqrt(eta)*u_mod)*u(a)*u(b)/u_mod
  !         if (a==b) then
  !            value=value+sqrt(eta)*h_(sqrt(eta)*u_mod) 
  !         end if
  !      end do  
  !   end do
  !end do
  !call stop_clock( 'real' )
  end subroutine add_local_uno

  subroutine I_due_value(y,x,flag)
  use cell_base, only:tpiba
  use gvect, only: ngm,gstart,g,gl,gg
  use mp,        only :mp_sum
  use wvfct,      only :npw 
  use mp_pools,   only : intra_pool_comm
  implicit none
  integer, intent(in) ::flag 
  real(DP), intent(out) ::y
  real(DP),intent(in)  ::x(3)
  integer ::igm!,ng_max
  real(DP) ::scalar
!
  call start_clock( 'rec' ) 
!  ng_max=10
  y=0.d0
  do igm=gstart,npw
     if   ( (gg(igm)*tpiba*tpiba)/(4.d0*eta) > 20.d0) exit
     y=y+2.d0*(I_due_g(igm)*cos(DOT_PRODUCT(g(1:3,igm), x(1:3)) *tpiba) )
!     if	(gl(igm)>ng_max) exit
  end do
  if (gstart==2) y=y+I_due_g(1)
  call mp_sum(y, intra_pool_comm)
  if (flag==1) then
     call add_local_due(y,x)
  end if
  call stop_clock( 'rec' )
  end subroutine i_due_value

  subroutine add_local_due(value,pos)
  use cell_base, only:at,alat
  use mp, only : mp_sum
  use mp_world,      ONLY : nproc,mpime
  use mp_pools,      ONLY : intra_pool_comm
  implicit none
  real(DP), intent(in) ::pos(3)
  real(DP), intent(inout) ::value
  real(DP) ::modul,n(3),erf_value
  integer  ::n_x,n_y,n_z
  real(kind=DP), external :: qe_erfc
  integer   :: l_blk,nbegin,nend
  real(DP)  :: value1
  
!
!!!!!!!!!!!
  call start_clock( 'real' )
  l_blk= (2*n_max+1)/nproc
  if(l_blk*nproc <  (2*n_max+1) ) l_blk = l_blk+1
  nbegin=mpime*l_blk-n_max
  nend=nbegin+l_blk-1
  if(nend > n_max) nend = n_max
  value1=0.d0
  do n_x=nbegin,nend
     do n_y=-n_max,n_max      !secondo ciclo su n
        do n_z=-n_max,n_max   !terzo ciclo   su n
           n(1:3) = n_x * at(1:3,1)*alat + n_y * at(1:3,2)*alat + n_z * at(1:3,3)*alat
           modul=modulus(pos(1:3)-n(1:3))
           erf_value=qe_erfc(sqrt(eta)*modul)
           value1=value1+ erf_value/modul
        end do
     end do
  end do
  call mp_sum(value1, intra_pool_comm)
  value=value+value1
  call stop_clock( 'real' )
!!!!!!!!!!
!  call start_clock( 'real' )
!  do n_x=-n_max,n_max         !primo ciclo   su n
!     do n_y=-n_max,n_max      !secondo ciclo su n
!        do n_z=-n_max,n_max   !terzo ciclo   su n                 
!           n(1:3) = n_x * at(1:3,1)*alat + n_y * at(1:3,2)*alat + n_z * at(1:3,3)*alat
!           modul=modulus(pos(1:3)-n(1:3)) 
!           erf_value=qe_erfc(sqrt(eta)*modul)
!           value=value+ erf_value/modul                      
!        end do  
!     end do
!  end do
!  call stop_clock( 'real' )

  end subroutine add_local_due

!!!!!!!!!!!!!!!!!!!
END MODULE zero_mod

