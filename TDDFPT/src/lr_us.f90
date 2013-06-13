!!----------------------------------------------------------------------------
module lr_us
!----------------------------------------------------------------------------
! Created by Xiaochuan Ge (May, 2013)
!-----------------------------------------------------------------------
contains

  subroutine lr_apply_s(vect,svect)
    !-------------------------------------------------------------------------------
    ! Created by X.Ge in May. 2013
    !-------------------------------------------------------------------------------
    ! This routine apply S to vect (input,unchanged at the exit), svect is output.

    use kinds,          only : dp
    use uspp,           only : okvan, vkb
    use wvfct,          only : npwx, npw, nbnd
    use klist,          only : nks,xk
    USE control_flags,  ONLY : gamma_only
    use becmod,         only : becp, calbec
    USE realus,               ONLY : real_space, fft_orbital_gamma, &
                                    bfft_orbital_gamma, calbec_rs_gamma, &
                                    v_loc_psir, igk_k,npw_k, real_space_debug

    implicit none
    complex(dp) :: vect(npwx,nbnd,nks),svect(npwx,nbnd,nks),vect_temp(npwx,nbnd,nks)
    integer :: ibnd,ik
    
    vect_temp(:,:,:)=vect(:,:,:)
    if( .not. okvan ) then  ! Not USPP
      svect(:,:,:)=vect(:,:,:)
    else ! USPP
      if(gamma_only) then
        if (real_space_debug>5) THEN ! real space & nkb > 0
          do ibnd=1, nbnd, 2
            call fft_orbital_gamma(vect(:,:,1),ibnd,nbnd)
            call calbec_rs_gamma(ibnd,nbnd,becp%r)
            call bfft_orbital_gamma(svect(:,:,1),ibnd,nbnd)
          enddo
        else
          CALL calbec(npw,vkb,vect(:,:,1),becp)
          CALL s_psi(npwx,npw,nbnd,vect(1,1,1),svect(1,1,1))
        endif
      else ! Generalised K points algorithm
        do ik = 1, nks
          CALL init_us_2(npw_k(ik),igk_k(1,ik),xk(1,ik),vkb)
          CALL calbec(npw_k(ik),vkb,vect(:,:,ik),becp)
          CALL s_psi(npwx,npw_k(ik),nbnd,vect(1,1,ik),svect(1,1,ik))
        enddo
      endif
    endif
   
    return
  end subroutine lr_apply_s
  !-------------------------------------------------------------------------------

  function lr_dot_us(vect1,vect2)
    !-------------------------------------------------------------------------------
    ! Created by X.Ge in May. 2013
    !-------------------------------------------------------------------------------
    ! This routine calculates < vect1 | S | vect2 >
    use kinds,          only : dp
    use wvfct,          only : npwx, npw, nbnd
    use klist,          only : nks

    implicit none
    complex(kind=dp),external :: lr_dot
    complex(dp) :: vect1(npwx,nbnd,nks),vect2(npwx,nbnd,nks),svect1(npwx,nbnd,nks)
    complex(dp) :: lr_dot_us

    call lr_apply_s(vect1(1,1,1),svect1(1,1,1))
    lr_dot_us=lr_dot(svect1,vect2)
    
    return
  end function lr_dot_us
  !-------------------------------------------------------------------------------

END MODULE lr_us
