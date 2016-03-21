subroutine qpcorrections(wcstate)
! this subroutine reads the qp corrections and gives the correct weight to each line of the 
! wcstate vector

!USE qpe,         ONLY: qpc,qpcbarc,qpcbarv
USE qpe_exc         
USE bse_wannier, ONLY: qpe_imin, qpe_imax, num_nbndv,scissor
USE kinds,        ONLY : DP
USE bse_basic_structures

implicit none
type(c_state)   :: wcstate
real(kind=DP), allocatable :: qpcw(:)


call build_qpc(qpc)
allocate(qpcw(wcstate%numb_c))
qpcw=0.d0

if (qpe_imin <= num_nbndv(1)) then
   qpcbarv=qpc(qpe_imin)
   qpc(1:qpe_imin)=qpcbarv
else
!case only conduction bands corrections are computed, valence shifts rigidly 
   qpc(1:qpe_imin)=-scissor
endif

if (qpe_imax > num_nbndv(1)) then
   qpcbarc=qpc(qpe_imax)
   qpcw(1:qpe_imax-num_nbndv(1))=qpc(num_nbndv(1)+1:qpe_imax)-qpcbarc
   
   call c_times_cstate(qpcw,wcstate,wcstate)
else
!case only valence bands corrections are computed, conduction shifts rigidly 
   qpcbarc=scissor
endif




return
end subroutine
