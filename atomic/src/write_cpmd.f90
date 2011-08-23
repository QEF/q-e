!
! Copyright (C) 2005 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
!  Translates from usual dft acronyms to appropriate indexes for CPMD
!-----------------------------------------------------------------------
subroutine which_cpmd_dft &
     (dft,mfxcx, mfxcc, mgcx, mgcc)
  !-----------------------------------------------------------------------
  !
  use funct, only : get_iexch, get_icorr, get_igcx, get_igcc, set_dft_from_name
  implicit none
  character(len=*), intent(IN) :: dft
  integer , intent(OUT) :: mfxcx, mfxcc, mgcx, mgcc

  call set_dft_from_name(dft)
  mfxcx = get_iexch()
  mfxcc = get_icorr()
  mgcx = get_igcx()
  mgcc = get_igcc()

! in CPMD PW91 and LYP are swapped.
  if (mgcc.eq.3) then
     mgcc=2
  else if (mgcc.eq.2) then
     mgcc=3
  end if

  return
end subroutine which_cpmd_dft
!
!-----------------------------------------------------------------------
subroutine write_cpmd &
     (iunps,zed,xmin,dx,mesh,ndm,r,r2, &
     dft,lmax,lloc,zval,nlc,nnl,cc,alpc,alc,alps,nlcc, &
     rhoc,vnl,phis,vpsloc,els,lls,ocs,rcuts,etots,nwfs)
  !-----------------------------------------------------------------------
  !
  use kinds, only : DP
  use constants, only : fpi, e2
  implicit none
  integer :: iunps, ndm, mesh, nwfps, lmin,lmax,lloc,nlc,nnl,nwfs,lls(nwfs)
  real(DP) :: zed, zval, xmin,dx, cc(2),alpc(2),alc(6,0:3), &
       alps(3,0:3), phis(ndm,nwfs), ocs(nwfs), rcuts(nwfs), &
       r(ndm), r2(ndm), vnl(ndm,0:3), rhoc(ndm), etots
  character(len=*) :: dft
  !
  real(DP) :: alfa_core=0.0_dp, vpsloc(ndm)
  logical nlcc, bhstype, numeric
  character(len=70) title_pseudo
  character(len=2), external :: atom_name
  character(len=2) :: els(nwfs)
  integer :: ios, i, l, k, n, ir, nb, mfxcx, mfxcc, mgcx, mgcc, lls_table(nwfs)
  !
  !
  nlc=0
  nnl=0
  bhstype=.false.

  call which_cpmd_dft(dft, mfxcx, mfxcc, mgcx, mgcc)

  do l=0,lmax
     do nb = 1, nwfs
        if (lls(nb).eq.l) then
           lls_table(l+1)=nb
        endif
     enddo
  enddo

  if (nlcc) then
     write(title_pseudo,"(2a3,'cc',f4.2,4(1x,a2,f4.2,' Rc=',f4.2))")  &
          'MT', atom_name(nint(zed)), alfa_core,                      &
          (els(lls_table(n)),ocs(lls_table(n)),rcuts(lls_table(n)),n=1,lmax+1)
  else
     write(title_pseudo,"(2a3,4(1x,a2,f4.2,' Rc=',f4.2))")  &
          'MT', atom_name(nint(zed)),                       &
          (els(lls_table(n)),ocs(lls_table(n)),rcuts(lls_table(n)),n=1,lmax+1)
  endif

  write(iunps, "('&ATOM')", err=300, iostat=ios)
  write(iunps, "(' Z =',i2)", err=300, iostat=ios) nint(zed)
  write(iunps, "(' ZV=',i2)", err=300, iostat=ios) nint(zval)
               
  write(iunps, "(' XC=',4i1,8x,'.666667')", err=300, iostat=ios) mfxcx,mfxcc,mgcx,mgcc
  write(iunps, "(' TYPE=NORMCONSERVING NUMERIC')", err=300, iostat=ios)
  write(iunps, "('&END')", err=300, iostat=ios)
  write(iunps, "('&INFO')", err=300, iostat=ios)
  write(iunps, '(1x,a)', err=300, iostat=ios) title_pseudo
  write(iunps, "('&END')", err=300, iostat=ios)

  write(iunps, "('&POTENTIAL')", err=300, iostat=ios)
  write(iunps, '(i6,f15.8)', err=300, iostat=ios) mesh,exp(dx)
  do i=1,mesh
     write(iunps, '(5e18.8)', err=300, iostat=ios) r(i),((vnl(i,l)+vpsloc(i))/e2,l=0,lmax)
  end do
  write(iunps,"('&END')", err=300, iostat=ios)

  write(iunps,"('&WAVEFUNCTION')", err=300, iostat=ios)
  write(iunps,'(i6,f15.8)', err=300, iostat=ios) mesh,exp(dx)
  do i=1,mesh
     write(iunps,'(5e18.8)', err=300, iostat=ios) r(i),(phis(i,lls_table(n)),n=1,lmax+1)
  end do
  write(iunps,"('&END')", err=300, iostat=ios)
  if(nlcc) then
     write(iunps,"('&NLCC')", err=300, iostat=ios)
     write(iunps,"('     NUMERIC')", err=300, iostat=ios)
     write(iunps,'(i4)', err=300, iostat=ios) mesh
     write(iunps,'(2e16.8)', err=300, iostat=ios) (r(i), rhoc(i), i=1,mesh)         
!         write(iunps,'(2e16.8)') (r(i), rho_core(i)/fpi, i=1,mesh)       
  
     write(iunps,"('&END')", err=300, iostat=ios)
  end if
  return

300 call errore('write_cpmd','writing pseudo file',abs(ios))

end subroutine write_cpmd
