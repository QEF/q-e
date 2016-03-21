MODULE qpe_exc
!this modules contains variables and subroutines related to the use of QP energies
!in the BSE kernel beyond the scissor operator
 
USE kinds, ONLY : DP
USE io_global, ONLY : ionode,ionode_id
USE mp, ONLY : mp_bcast
USE mp_world, ONLY : world_comm

REAL(kind=DP), pointer :: qpc(:)  ! vector containing QPC
REAL(kind=DP)          :: qpcbarc ! average qpc to be applied for higher 
                                  ! (band index above qpc_imax) energy states
REAL(kind=DP)          :: qpcbarv ! average qpc to be applied for lower 
                                  ! (band index below qpc_imin) energy states

CONTAINS

    SUBROUTINE build_qpc(qpc)
    !this subroutine reads bands.dat and builds the qp correction vector
    USE bse_wannier, ONLY : qpe_imin,qpe_imax
    USE kinds, ONLY : DP
    USE io_files, ONLY : tmp_dir,prefix
    USE wvfct,    ONLY : nbnd
    USE constants,   ONLY: RYTOEV

    implicit none
     
    INTEGER, EXTERNAL :: find_free_unit
    integer :: ib, iun, idum
    logical :: debug
    real(kind=DP) :: qpc(qpe_imax)
    real(kind=DP)          :: rdum,edft,egw,dumm1,dumm2,dumm3,dumm4
        
    qpc(1:qpe_imax)=0.d0

    iun = find_free_unit()  
    if(ionode) then
       open(iun,file=trim(tmp_dir)//trim(prefix)//'-bands.dat', status='old', form='formatted')
       
       read(iun,*) idum
       read(iun,*) idum

       do ib=1,qpe_imin-1
          read(iun,*) idum,dumm1,dumm2,dumm3,dumm4
       enddo

       do ib=qpe_imin,qpe_imax
          read(iun,*) idum, edft,rdum,egw,rdum
          qpc(ib)=(egw-edft)/RYTOEV
       
       enddo
 
       close(iun)
    endif
    do ib=qpe_imin,qpe_imax
       call mp_bcast(qpc(ib), ionode_id, world_comm )
    enddo
    return
    END SUBROUTINE

end MODULE qpe_exc
