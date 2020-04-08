subroutine build_eemat(data_input,bd,eemat)
  USE input_simple_exc
  USE simple_objects
  USE derived_objects
  USE kinds, ONLY : DP
  USE io_global, ONLY : stdout, ionode
  USE constants, ONLY : rytoev
  USE mp, ONLY : mp_barrier, mp_sum
  USE mp_world,             ONLY : world_comm,mpime
  use io_files, only : prefix, tmp_dir
  USE constants, ONLY : rytoev

  implicit none

  TYPE(input_options) :: data_input
  TYPE(epe) :: element
  TYPE(bands) :: bd
  COMPLEX(KIND=DP),ALLOCATABLE::Etempc(:,:,:),Etempv(:,:,:),mate(:,:,:,:)
  COMPLEX(KIND=DP) :: eemat(bd%numv,bd%numc,bd%nk_loc,3)
  INTEGER :: a, i, j, v, c, k
  
  allocate(Etempv(bd%ntot_e,bd%numv,bd%nk_loc))!vettori E_{iv}^k
  allocate(Etempc(bd%ntot_e,bd%numc,bd%nk_loc))!vettori E_{ic}^k
  allocate(mate(bd%ntot_e,bd%numc,bd%nk_loc,3))!E_{jc}^k <e_i|p_a|e_j>

  call initialize_epe(element)
  write(stdout,*)'reading epe'
  call read_epe(data_input,element)
  write(stdout ,*)'EPE MATRIX DIMENSION', element%ntot_e
  

  !costruisco E_{jc}^k
  do k=1,bd%nk_loc
     do c=bd%numv+1,bd%num
        do j=1,bd%ntot_e
           Etempc(j,c-bd%numv,k) = bd%omat(j,c,k)
        end do
     end do
  end do

  write(stdout,*) 'E_{jc}^k'

  !prodotto <e_i|p_a|e_j>E_{jc}^k
  do a=1,3
     do k=1,bd%nk_loc
        call zgemm('N','N',bd%ntot_e,bd%numc,bd%ntot_e,(1.d0,0.d0),element%mat(:,:,a),bd%ntot_e,Etempc(:,:,k), &
             &bd%ntot_e,(0.d0,0.d0),mate(:,:,k,a),bd%ntot_e)
     end do
  end do
  write(stdout,*)'E_{jc}^k <e_i|p_a|e_j>'
  !costruisco E_{iv}^k
  do k=1,bd%nk_loc
     do v=1,bd%numv
        do i=1,bd%ntot_e
           Etempv(i,v,k) = bd%omat(i,v,k)
        end do
     end do
  end do
  write(stdout,*) 'E_{iv}^k'
  !prodotto E_{iv}^k^*  E_{jc}^k <e_i|p_a|e_j>
  do a=1,3
     do k=1,bd%nk_loc
        call zgemm('C','N',bd%numv,bd%numc,bd%ntot_e,(1.d0,0.d0),Etempv(1,1,k),bd%ntot_e,mate(1,1,k,a),bd%ntot_e,(0.d0,0.d0),&
             &eemat(1,1,k,a),bd%numv)
     end do
  end do
  write(stdout,*) ' E_{iv}^k^*  E_{jc}^k <e_i|p_a|e_j>'

  deallocate(Etempv)
  deallocate(Etempc)
  deallocate(mate)
  call deallocate_epe(element)

  
  return
end subroutine build_eemat
