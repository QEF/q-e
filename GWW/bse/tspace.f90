subroutine tspace_diago(vstate,vstate_r,fc)
! this subroutine solves teh BSE equation in transition space

USE bse_wannier, ONLY : num_nbndv
USE transitions, ONLY : ttab,itiv,itic,exch
USE wvfct,       ONLY : nbnd
USE io_global,   ONLY : stdout,ionode
USE kinds,       ONLY : DP
USE constants,   ONLY : RYTOEV
USE mp,          ONLY : mp_barrier
USE mp_world,             ONLY : world_comm
USE bse_basic_structures 
USE fft_custom_gwl 


implicit none
type(v_state) vstate
type(v_state_r) vstate_r
type(fft_cus) :: fc

integer                         :: nt,it ! number of transitions
real(kind=DP),    allocatable   :: eig(:)

! for zheev
integer                         :: info
real(kind=DP), allocatable   :: work(:)

call start_clock('tspace_diago')
! compute the number of transitions
nt=num_nbndv(1)*(nbnd-num_nbndv(1))
if(ionode) write(stdout,*) 'number of transitions nt:', nt

! allocate and build the transition table (c,v)-->it
allocate(ttab(nbnd,num_nbndv(1)))
allocate(itiv(nt))
allocate(itic(nt))

call build_ttab

! allocate and build the excitonic Hamiltonian in transition space
allocate(exch(nt,nt))
exch(1:nt,1:nt)=0.d0

call mp_barrier(world_comm)
call build_exch(vstate,vstate_r,fc)

call mp_barrier(world_comm)
! diagonalize (not parallel so only one proc will do it)

allocate(eig(nt))
if(ionode) then
   allocate(work(3*nt-1))
   call dsyev('V', 'U', nt, exch, nt, eig, work, 3*nt-1 , info)
   eig(1:nt)=eig(1:nt)*RYTOEV 
   deallocate(work)
endif

call mp_barrier(world_comm)

if(ionode) then
   do it=1,nt 
      write(stdout,*) 'Eigenvalue number', it, eig(it)
   enddo
endif

deallocate(ttab,itic,itiv,exch,eig)
call stop_clock('tspace_diago')
return
end subroutine

!--------------------!
subroutine build_exch(vstate,vstate_r,fc)
!this subroutine builds the excitonic Hamiltonian in transition space

USE bse_wannier,           ONLY : num_nbndv
USE wvfct,                 ONLY : nbnd, npwx,npw,et
USE lsda_mod,              ONLY : nspin
USE wavefunctions_module,  ONLY : evc
USE io_files,              ONLY : prefix, iunwfc, nwordwfc
USE transitions,           ONLY : ttab,exch
USE io_global,             ONLY : stdout,ionode
USE bse_basic_structures 
USE fft_custom_gwl 
USE exciton 


implicit none
type(v_state),intent(in)    :: vstate
type(v_state_r), intent(in) :: vstate_r
type(fft_cus), intent(in)   :: fc

!internal
type(exc) :: phic
type(exc) :: phicp
type(exc) :: hphic
integer   :: iv,ivp,ic,icp,is

call start_clock('build_exch')

call initialize_exc(phic)
phic%label=1
phic%npw=npw
phic%numb_v=num_nbndv(1)
allocate(phic%a(phic%npw,phic%numb_v))

call initialize_exc(hphic)
hphic%label=1
hphic%npw=npw
hphic%numb_v=num_nbndv(1)
allocate(hphic%a(hphic%npw,hphic%numb_v))

call initialize_exc(phicp)
phicp%label=1
phicp%npw=npw
phicp%numb_v=num_nbndv(1)
allocate(phicp%a(phicp%npw,phicp%numb_v))

allocate( evc( npwx, nbnd ) )

!read wavefunctions
do is=1,nspin
   call davcio(evc,2*nwordwfc,iunwfc,is,-1)
enddo

if(ionode) write(stdout,*) 'wfns read from disk'

do iv=1,num_nbndv(1)
   do ic=num_nbndv(1)+1,nbnd
      
      if(ionode) write(stdout,*) 'iv ic', iv, ic
      phic%a(1:phic%npw,1:phic%numb_v)=dcmplx(0.d0,0.d0)
      phic%a(1:phic%npw,iv)=evc(1:npw,ic)
      call normalize_exc(phic)
      !apply tyhe excitonic Hamiltonian to this ic,iv state
      call exc_h_a(phic,hphic,vstate,vstate_r,fc) 
     
      do ivp=iv,num_nbndv(1)
          do icp= ic,nbnd
             if(ionode) write(stdout,*) 'ivp icp', ivp, icp
             phicp%a(1:phicp%npw,1:phicp%numb_v)=dcmplx(0.d0,0.d0)
             phicp%a(1:phicp%npw,ivp)=evc(1:npw,icp)
             call normalize_exc(phicp)
             call sproduct_exc(phicp,hphic,exch(ttab(ic,iv),ttab(icp,ivp)))
             if(ionode) write(stdout,*) 'Exc. Hamiltonian built for transition:', ttab(ic,iv),ttab(icp,ivp) 
          enddo !icp
      enddo !ivp
   enddo !ic
enddo !iv

call free_memory_exc_a(phic)
deallocate(evc)

call stop_clock('build_exch')
return
endsubroutine

!--------------------!
subroutine build_ttab
!this subroutine builds the transition table

USE wvfct,    ONLY : nbnd
USE bse_wannier, ONLY:num_nbndv
use transitions, ONLY:ttab,itiv,itic
USE io_global, ONLY : stdout,ionode

implicit none
integer       :: it,iv,ic

it=1

do iv=1,num_nbndv(1)
   do ic=num_nbndv(1)+1,nbnd
      ttab(ic,iv)=it
      itiv(it)=iv
      itic(it)=ic
      it=it+1
   enddo
enddo

it=it-1

if(ionode) write(stdout,*) 'ttab built, number of transitions found:', it
if(ionode) write(stdout,*) 'total number of bands', nbnd
if(ionode) write(stdout,*) 'number of valence bands', num_nbndv(1)
return
end subroutine
!--------------------!
