!
! Copyright (C) 2001-2015 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------
module lr_dav_debug
!----------------------------------------------------------------------------
! Created by Xiaochuan Ge (Aprile, 2013)
!-----------------------------------------------------------------------
contains

  ! subroutine for debug
  subroutine check_orth()
    !-------------------------------------------------------------------------------
    ! Created by X.Ge in Jan. 2013
    !-------------------------------------------------------------------------------
    use lr_dav_variables
    use io_global,   only : stdout
    use lr_us,   only : lr_dot_us
    use uspp,    only : okvan
    
    implicit none
    complex(kind=dp),external :: lr_dot
    integer :: ia, ib
    real(dp) :: inner_err,temp
    
    inner_err=0.0d0
    do ia = 1, num_basis
      do ib = 1, num_basis 
        if( poor_of_ram .or. .not. okvan) then
          inner_matrix(ia,ib)=dble(lr_dot_us(vec_b(:,:,1,ia),vec_b(:,:,1,ib)))
        else
          inner_matrix(ia,ib)=dble(lr_dot(svec_b(:,:,1,ia),vec_b(:,:,1,ib)))
        endif
        if(ia==ib) then
          temp=(inner_matrix(ia,ib)-1)**2
        else
          temp=inner_matrix(ia,ib)**2
        endif
        inner_err=inner_err+temp
        if(temp .gt. 10E-10) &
          &write(stdout,*) "Warning, the inner product between ",ia," and ",ib," is : ",temp
      enddo
    enddo
    
    inner_err=inner_err/(num_basis*num_basis)
    write(stdout,'(/5x,"The error of the orthonalization of the basis is:",&
            & 5x,E20.12)') inner_err
    call check("inner_matrix")
    return
  end subroutine check_orth
  !-------------------------------------------------------------------------------

  subroutine check(flag_check)
    !-------------------------------------------------------------------------------
    ! Created by X.Ge in Jan. 2013
    !-------------------------------------------------------------------------------
    ! For debuging
    use lr_dav_variables
    use kinds,     only : dp
    use io_global, only : stdout
    USE wvfct,         ONLY : nbnd,npwx
    USE klist,             ONLY : nks
    use lr_variables,    only : evc0, sevc0

    implicit none
    integer :: ia,ib
    character(len=*) :: flag_check
   
    if(.not. dav_debug) return

    write(stdout,'(/5x,"You are checking some information for debugging")')
    if(flag_check == "M_C") then
      write(stdout, '(7x,"Matrix C is:")')
      do ia=1,num_basis
        write(stdout,'(7x,10F15.8)') (dble(M_C(ia,ib)),ib=1,min(10,num_basis))
      enddo
    endif

    if(flag_check == "M_D") then
      write(*, '(7x,"Matrix D is:")')
      do ia=1,num_basis
        write(stdout,'(7x,10F15.8)') (dble(M_D(ia,ib)),ib=1,min(10,num_basis))
      enddo
    endif

    if(flag_check == "M") then
      write(*, '(7x,"Matrix DC is:")')
      do ia=1,num_basis
        write(stdout,'(7x,10F15.8)') (dble(M(ia,ib)),ib=1,min(10,num_basis))
      enddo
    endif

    if(flag_check == "right_M") then
      write(*, '(7x,"Matrix right_M is:")')
      do ia=1,num_basis
        write(stdout,'(7x,10F15.8)') (dble(right_M(ia,ib)),ib=1,min(10,num_basis))
      enddo
    endif

    if(flag_check == "left_M") then
      write(*, '(7x,"Matrix left_M is:")')
      do ia=1,num_basis
        write(stdout,'(7x,10F15.8)') (dble(left_M(ia,ib)),ib=1,min(10,num_basis))
      enddo
    endif

    if(flag_check == "inner_matrix") then
      write(*, '(7x,"Inner_matrix is:")')
      do ia=1,num_basis
        write(stdout,'(7x,10F15.8)') (dble(inner_matrix(ia,ib)),ib=1,min(10,num_basis))
      enddo
    endif

    return
  end subroutine check
  !-------------------------------------------------------------------------------

  subroutine check_overlap(vec)
    !-------------------------------------------------------------------------------
    ! Created by X.Ge in Apr. 2013
    !-------------------------------------------------------------------------------
    ! To see if the vector is othogonal to occupied space
    
    use lr_variables,         only : evc0
    use wvfct,                only : nbnd,npwx
    use io_global,            only : stdout
    use kinds,                 only : dp
    use lr_us,     only : lr_dot_us

    implicit none
    complex(dp) :: vec(npwx,nbnd),occupy(npwx,nbnd), overlap
    integer :: i,j

    overlap=0.0d0
    do i = 1, nbnd
      do j = 1, nbnd
        occupy(:,j)=evc0(:,i,1)
      enddo
      overlap=overlap+lr_dot_us(vec(:,:),occupy(:,:))
    enddo
    write(stdout,'("!!!! the tot overlap with the occupied space is:",5x,E20.12)') dble(overlap)

    return
  end subroutine check_overlap
  !-------------------------------------------------------------------------------

  subroutine check_overlap_basis(vec)
    !-------------------------------------------------------------------------------
    ! Created by X.Ge in Apr. 2013
    !-------------------------------------------------------------------------------
    ! Check the overlap between residue and basis
    use lr_variables,         only : evc0
    use wvfct,                only : nbnd,npwx
    use kinds,                 only : dp
    use io_global,            only : stdout
    use lr_dav_variables,      only : vec_b,num_basis

    implicit none

    complex(kind=dp),external :: lr_dot
    complex(dp) :: vec(npwx,nbnd),overlap
    integer :: i

    overlap=0.0d0
    do i = 1, num_basis
      overlap = overlap + lr_dot(vec(:,:),vec_b(:,:,1,i))**2
    enddo

    write(stdout,'("!!!! the tot overlap of the residue with the basis space is:",5x,E20.12)') dble(overlap)
   
    return
  end subroutine check_overlap_basis
  !-------------------------------------------------------------------------------

  subroutine check_revc0()
    !-------------------------------------------------------------------------------
    ! Created by X.Ge in Apr. 2013
    !-------------------------------------------------------------------------------
    ! Due to the bug of virt_read, this is to check if revc0 is correct
    use lr_variables,    only : evc0, revc0
    use wvfct,           only : nbnd, npwx
    use kinds,           only : dp
    use fft_base,             only : dffts
    use mp_global,            only : intra_bgrp_comm
    use mp_world,             only : world_comm
    use mp,                   only : mp_sum, mp_barrier
    use lr_dav_variables
    USE cell_base,              ONLY : omega
    USE wavefunctions_module, ONLY : psic
    USE realus,              ONLY : invfft_orbital_gamma, fwfft_orbital_gamma
      USE gvect,                ONLY : gstart

    implicit none
    integer :: i,tot_nnr
    real(dp) :: norm,banda(dffts%nnr),bandb(dffts%nnr)
    real(kind=dp), external    :: ddot
    complex(kind=dp),external :: lr_dot   
    complex(dp) :: wfck(npwx,1)

    ALLOCATE( psic(dffts%nnr) )
    
    tot_nnr=dffts%nr1x*dffts%nr2x*dffts%nr3x
    
    do i = 1, nbnd,2
      banda(:)=dble(revc0(:,i,1))     
      bandb(:)=aimag(revc0(:,i,1))     
      wfck(:,1)=evc0(:,i,1)
      call invfft_orbital_gamma(wfck(:,:),1,1)  ! FFT: v  -> psic
      norm=DDOT(dffts%nnr,psic(:),2,banda,1)/tot_nnr
#if defined(__MPI)
      call mp_barrier( world_comm )
      call mp_sum(norm,intra_bgrp_comm)
#endif
      print *, norm

      wfck(:,1)=evc0(:,i+1,1)
      call invfft_orbital_gamma(wfck(:,:),1,1)  ! FFT: v  -> psic
      norm=DDOT(dffts%nnr,psic(:),2,bandb,1)/tot_nnr
#if defined(__MPI)
      call mp_barrier( world_comm )
      call mp_sum(norm,intra_bgrp_comm)
#endif
      print *, norm
    enddo
!call mp_stop(100)
    return
  end subroutine check_revc0
  !-------------------------------------------------------------------------------

  subroutine check_hermitian()
    !-------------------------------------------------------------------------------
    ! Created by X.Ge in Apr. 2013
    !-------------------------------------------------------------------------------
    use lr_dav_variables
    use lr_variables, only : evc0,sevc0
    use kinds,     only : dp
    use lr_us,  only : lr_dot_us

    implicit none
    integer :: ibr,ibl,ibr0,ibl0
    real(dp) :: diff,diffmax

    do ibr = 1, num_basis
      call lr_apply_liouvillian(vec_b(:,:,:,ibr),vecwork(:,:,:),.true.)
      do ibl = 1, num_basis
        M_C(ibl,ibr)=lr_dot_us(vec_b(1,1,1,ibl),vecwork(1,1,1))
      enddo
    enddo

    diffmax=0.0d0
    ibr0=1;ibl0=1
    do ibr = 1, num_basis-1
      do ibl = ibr+1, num_basis
        diff = abs(dble(M_C(ibl,ibr))-dble(M_C(ibr,ibl)))
        print *, diff,dble(M_C(ibl,ibr)),dble(M_C(ibr,ibl))
        if( diff .gt. diffmax ) then
          diffmax=diff
          ibr0=ibr;ibl0=ibl
        endif
      enddo
    enddo

    print *, "Max|C(i,j)-C(j,i)|=",diffmax
    print *, "i,j=",ibl0,ibr0
    print *, "C(i,j);C(j,i)",dble(M_C(ibl0,ibr0)),dble(M_C(ibr0,ibl0))

  end subroutine check_hermitian
  !-------------------------------------------------------------------------------

END MODULE lr_dav_debug
