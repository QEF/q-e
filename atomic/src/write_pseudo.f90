!
! Copyright (C) 2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine write_pseudo &
     (iunps,zed,xmin,dx,mesh,ndm,r,r2, &
     dft,lmax,lloc,zval,nlc,nnl,cc,alpc,alc,alps,nlcc, &
     rhoc,vnl,phis,vpsloc,els,lls,ocs,etots,nwfs)
  !-----------------------------------------------------------------------
  !
  use kinds, only : DP
  use constants, only : fpi
  use funct , only : dft_is_nonlocc 
  implicit none
  integer :: ndm, mesh, nwfps, lmin,lmax,lloc,nlc,nnl,nwfs,lls(nwfs)
  real(DP) :: zed, zval, xmin,dx, cc(2),alpc(2),alc(6,0:3), &
       alps(3,0:3), phis(ndm,nwfs), ocs(nwfs), &
       r(ndm), r2(ndm), vnl(ndm,0:3), rhoc(ndm), etots
  integer ios, mdum, i, l, k, n, ir, iunps, nb, ldum
  real(DP) :: zdum, vnloc,a_core, b_core,  &
       alfa_core, dum, xdum, dxdum, rdum, vpsloc(ndm)
  logical nlcc, bhstype, numeric
  character(len=70) title_pseudo
  character(len=2), external :: atom_name
  character(len=2) :: els(nwfs)
  character(len=*) dft
  logical :: non_locc
  !
  !
  nlc=0
  nnl=0
  bhstype=.false.
  non_locc = dft_is_nonlocc()

  if (dft.eq.'PW') then
     write( iunps, '(a)', err=300, iostat=ios ) 'slater-pz-ggx-ggc'
  else
     if (non_locc) then
        CALL errore('write_pseudo','non-local functionals not implemented yet', 1) 
     else
         write( iunps, '(a)', err=300, iostat=ios ) dft(1:20)
     endif

  endif

  write ( iunps, '("''",a2,"''",f8.4,3i5,l4,i5,l4,e17.9)', &
       err=300, iostat=ios )  atom_name(nint(zed)), &
       zval, lmax, nlc, nnl, nlcc, &
       lloc, bhstype, etots

  !
  !   In numeric pseudopotentials both nlc and nnl are zero.
  !
  numeric = nlc.le.0 .and. nnl.le.0

  if (.not.numeric) then
     write( iunps, *, err=300, iostat=ios ) &
          ( alpc(i), i=1, 2 ), ( cc(i), i=1,2 )
     do l = 0, lmax
        write ( iunps, *, err=300, iostat=ios ) &
             ( alps(i,l),i=1,3 ), (alc(i,l),i=1,6)
     enddo
     if (nlcc) then
        write( iunps, *, err=300, iostat=ios ) a_core, &
             b_core, alfa_core
     endif
  endif

  write( iunps, '(3f15.10,2i5)', err=300, iostat=ios ) &
       zed, xmin, dx, mesh, nwfs

  if (numeric) then
     !
     !      pseudopotentials in numeric form
     !
     do l = 0, lmax
        write( iunps, "(' Pseudopot. l=',i1)") l
        write( iunps, '(4e19.11)', err=300, iostat=ios ) &
             (vnl(ir,l)+vpsloc(ir),ir=1,mesh)
     enddo
     if (lloc.eq.-1) then
        write( iunps, "(' Local PP')", err=300, iostat=ios )
        write( iunps, '(4e19.11)', err=300, iostat=ios ) &
             (vpsloc(ir),ir=1,mesh)
     endif
     if(nlcc) then
        write( iunps, '(4e19.11)', err=300, iostat=ios ) &
             ( rhoc(ir)/r2(ir)/fpi, ir=1,mesh )
     endif
  endif

  do l=0,lmax
     do nb = 1, nwfs
        if (lls(nb).eq.l) then
           write( iunps, "(' Wavefunction ',a2)", err=300, iostat=ios ) els(nb)
           write( iunps, *, err=300, iostat=ios ) lls(nb), &
                (abs(ocs(nb))+ ocs(nb))*0.5_dp 
           write( iunps, '(4e19.11)', err=300, iostat=ios ) &
                (phis(ir,nb),ir=1,mesh)
        endif
     enddo
  enddo
300 call errore('write_pseudo','writing pseudo file',abs(ios))

  return
end subroutine write_pseudo
