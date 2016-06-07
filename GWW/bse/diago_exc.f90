subroutine diago_exc(a,v,cstate,wcstate)
! this subroutine applies the diagonal part of the excitonic Hamiltonian to the excitonic
! wavefunction vector (a%a)

USE kinds,            ONLY : DP
USE exciton
use bse_basic_structures
USE wvfct,            ONLY : g2kin,npwx
USE noncollin_module, ONLY : npol
USE uspp,             ONLY : vkb,nkb
USE becmod,           ONLY : becp,allocate_bec_type,deallocate_bec_type
USE g_psi_mod,        ONLY : h_diag, s_diag
USE klist,            ONLY : xk, igk_k
USE gvect
USE cell_base,        ONLY: tpiba,tpiba2
USE constants,        ONLY: RYTOEV 
use io_global,        ONLY : stdout, ionode 
use bse_wannier,      ONLY : scissor,l_scissor,qpe_imin,qpe_imax
use qpe_exc


implicit none

type(exc) :: a
type(v_state) :: v
type(c_state) :: cstate
type(c_state) :: wcstate


type(exc) :: a1,a2
COMPLEX(kind=DP), ALLOCATABLE :: psi_1(:,:)
COMPLEX(kind=DP), ALLOCATABLE :: u_0(:,:)



logical :: debug
real(kind=dp) :: prod
real(kind=dp), allocatable :: vb_en(:)
integer :: is

call start_clock('diago_exc')
debug=.false.

allocate(psi_1(a%npw,a%numb_v))
allocate(u_0(a%npw,a%numb_v))

ALLOCATE( h_diag( npwx,npol ) )
ALLOCATE( s_diag( npwx,npol ) )

ALLOCATE(vb_en(a%numb_v))

!just copy a in a temporary variable to apply the different part of the diago
!Hamiltonian 

call initialize_exc(a1)
call initialize_exc(a2)

allocate(a1%a(a%npw,a%numb_v))

a1%npw=a%npw
a1%numb_v=a%numb_v
a1%label=20

allocate(a2%a(a%npw,a%numb_v))

a2%npw=a%npw
a2%numb_v=a%numb_v
a2%label=30

a2%a(1:a2%npw,1:a2%numb_v)=a%a(1:a%npw,1:a%numb_v)

call allocate_bec_type ( nkb, a%numb_v, becp)

IF ( nkb > 0 )  CALL init_us_2( a%npw, igk_k(1,1), xk(1,1), vkb )
g2kin(1:a%npw) = ( (g(1,igk_k(1:a%npw,1)) )**2 + &
       ( g(2,igk_k(1:a%npw,1)) )**2 + &
       ( g(3,igk_k(1:a%npw,1)) )**2 ) * tpiba2


psi_1(1:a%npw,1:a%numb_v)=a%a(1:a%npw,1:a%numb_v)

    
!calculate H|\phi_i>
call h_psi( a%npw, a%npw, a%numb_v,psi_1(1,1), u_0 )
a1%a(1:a%npw,1:a%numb_v)=u_0(1:a%npw,1:a%numb_v)


!project into the conduction manifold
do is = 1,v%nspin
  call pc_operator_exc(a1,v,is)
enddo

!check if everything is ok, the 'scalar' product of a1%a with a%a should be
!greater than e_lumo

if (debug) then
  call sproduct_exc(a,a1,prod)
  prod=prod*RYTOEV
  if(ionode) write(stdout,*) 'exc_diago, prod (eV)=',prod
  if(ionode) write(stdout,*) 'prod should be greater than LUMO level'
  FLUSH ( stdout )
end if

if(.not.l_scissor) then
   if (qpe_imax>a%numb_v) then
      do is=1,a%numb_v
         vb_en(is)= qpcbarc
      enddo
      call c_times_exc(a2,vb_en)
      a1%a(1:a%npw,1:a%numb_v)=a1%a(1:a%npw,1:a%numb_v)+a2%a(1:a%npw,1:a%numb_v)
      a2%a(1:a2%npw,1:a2%numb_v)=a%a(1:a%npw,1:a%numb_v) 
      call poutcstate_exc(a2,a2,cstate,wcstate)
      a1%a(1:a%npw,1:a%numb_v)=a1%a(1:a%npw,1:a%numb_v)+a2%a(1:a%npw,1:a%numb_v)
   else 
      do is=1,a%numb_v
         vb_en(is)= qpcbarc
      enddo
      call c_times_exc(a2,vb_en)
      a1%a(1:a%npw,1:a%numb_v)=a1%a(1:a%npw,1:a%numb_v)+a2%a(1:a%npw,1:a%numb_v)
   endif
endif



!multiply each line of the excitonic wavefunction vector with the corresponding
!single particle valence state energy

if(l_scissor) then
   do is=1,a%numb_v
      vb_en(is)= v%esp(is,1)-scissor 
   enddo

   call c_times_exc(a,vb_en)
   if (debug) then
     do is=1,a%numb_v
        prod=vb_en(is)*RYTOEV
        if(ionode) write(stdout,*) 'exc_diago, band i  (eV)=',prod
     enddo
   end if
else !not scissor
   do is=1,a%numb_v
      vb_en(is)= v%esp(is,1)+qpc(is)
   enddo
   call c_times_exc(a,vb_en)
endif


! sum-up the two terms
a%a(1:a%npw,1:a%numb_v)=-a%a(1:a%npw,1:a%numb_v)+a1%a(1:a%npw,1:a%numb_v)

deallocate(psi_1)
deallocate(u_0)

deallocate(h_diag)
deallocate(s_diag)

deallocate(vb_en)

call deallocate_bec_type(becp)
call free_memory_exc_a(a1)
call free_memory_exc_a(a2)

call stop_clock('diago_exc')
return
end subroutine

