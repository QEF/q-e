program xctest
  USE mp, ONLY: mp_start, mp_end
  use kinds, only: dbl
  use funct
  implicit none
  integer, parameter :: nnr = 1000
  integer, parameter :: nspin = 2

  
  CALL mp_start()

  iexch=1
  icorr=3
  igcx=1
  igcc=3
  CALL test_xc( nnr, nspin )
  
  CALL mp_end()
end program

subroutine test_xc( nnr,nspin )
  use kinds, only: dbl
  use funct
  implicit none
  integer, intent(in) :: nnr, nspin
  real(dbl) :: rhor( nnr, nspin )
  real(dbl) :: grhor( nnr, 3, nspin )
  real(dbl) :: rhon( nnr, nspin )
  real(dbl) :: grhon( nnr, 3, nspin )
  real(dbl) :: exc, excn
  integer :: ir, is, ipol
  do is = 1, nspin
    do ir = 1, nnr
      rhor( ir, is ) = DBLE( is * ir ) / DBLE( nspin * nnr )
      grhor( ir, 1, is ) = 0.13 * rhor( ir, is )
      grhor( ir, 2, is ) = 0.43 * rhor( ir, is )
      grhor( ir, 3, is ) = 0.73 * rhor( ir, is )
    end do
  end do
  rhon = rhor
  grhon = grhor
  !
  ! original CP xc selection
  !
      if (iexch==1.and.icorr==1.and.igcx==0.and.igcc==0) then
         ! LDA (Perdew-Zunger)
         call expxc(nnr,nspin,rhor,exc)
      else if (iexch==1.and.icorr==4.and.igcx==2.and.igcc==2) then
         ! PW91
         call ggapwold(nnr,nspin,grhor,rhor,exc)
      else if (iexch==1.and.icorr==3.and.igcx==1.and.igcc==3) then
         ! BLYP
         call ggablyp4(nnr,nspin,grhor,rhor,exc)
      else if (iexch==1.and.icorr==4.and.igcx==3.and.igcc==4) then
         ! PBE
         call ggapbe(nnr,nspin,grhor,rhor,exc)
      else
         call errore('exc-cor','no such exch-corr',1)
      end if
  !
  ! Wrapper to PW xc selection
  !
  call exch_corr_cp(nnr,nspin,grhon,rhon,excn)
  !
  write(6,*) 'EXC = ', exc, excn
  write(6,*) 'VXC '
  do is = 1, nspin
    do ir = 1, nnr
      WRITE(6,100) ir,is,rhor( ir, is ),rhon( ir, is )
    end do
  end do
  do ipol = 1, 3
  write(6,*) 'IPOL ', ipol
  do is = 1, nspin
    do ir = 1, nnr
      WRITE(6,100) ir,is,grhor( ir, ipol, is ),grhon( ir, ipol, is )
    end do
  end do
  end do
100 FORMAT( I5, I2, 1X, E15.8, 1X, E15.8 )
end subroutine
