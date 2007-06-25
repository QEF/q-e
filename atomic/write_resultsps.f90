!
! Copyright (C) 2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------
subroutine write_resultsps 
  !--------------------------------------------------------------
  use constants, only : rytoev
  use ld1inc
  use funct, only: get_dft_name
  implicit none

  integer :: i, j, n, m, l, ios
  real(DP) :: work(ndm), int_0_inf_dr, ravg, sij
  character (len=20) :: dft_name
  !
  !
  dft_name = get_dft_name()
  write(6,110)
110 format (/,5x,22('-'),' Testing the pseudopotential ',22('-'),/)
  write(6,1150) title
  if(rel.eq.1) write(6,'(5x,''scalar relativistic calculation'')')
  if(rel.eq.2) write(6,'(5x,''dirac relativistic calculation'')')
1150 format(5x,a75)
  write(6,1250) zed, zval
1250 format(/5x,'atomic number is',f6.2,'   valence charge is',f6.2)
  write(6,2300) dft_name(1:len_trim(dft_name)),lsd,isic,latt,beta,tr2
2300 format(5x,'dft =',a,'   lsd =',i1,' sic =',i1,' latt =',i1, &
       '  beta=',f4.2,' tr2=',1pe7.1)
  write(6,1270) mesh,r(mesh),xmin,dx
1270 format(5x,'mesh =',i4,' r(mesh) =',f10.5,' xmin =',f6.2,' dx =',f8.5)
  if (rel.lt.2) then
     write(6,1000)
1000 format(/5x,'n l     nl             e AE (Ry) ',  &
          '       e PS (Ry)    De AE-PS (Ry) ')
     write(6,1100) &
          (nnts(n),llts(n),elts(n),iswts(n),octs(n), &
          enl(nstoaets(n)),enlts(n), &
          enl(nstoaets(n))-enlts(n),  n=1,nwfts)
     write(13,1100)  &
          (nnts(n),llts(n),elts(n),iswts(n),octs(n), &
          enl(nstoaets(n)),enlts(n),  &
          enl(nstoaets(n))-enlts(n),  n=1,nwfts)
1100 format(4x,2i2,5x,a2,i4,'(',f5.2,')',f15.5,f15.5,f15.5)
  else
     write(6,1001)
1001 format(/5x,'n l  j  nl             e AE (Ry)',  &
          '       e PS (Ry)    De AE-PS (Ry) ')
     write(6,1101) &
          (nnts(n),llts(n),jjts(n),elts(n),iswts(n),octs(n), &
          enl(nstoaets(n)),enlts(n), &
          enl(nstoaets(n))-enlts(n),  n=1,nwfts)
     write(13,1101)  &
          (nnts(n),llts(n),jjts(n),elts(n),iswts(n),octs(n), &
          enl(nstoaets(n)),enlts(n),  &
          enl(nstoaets(n))-enlts(n),  n=1,nwfts)
1101 format(4x,2i2,f4.1,1x,a2,i4,'(',f5.2,')',f15.5,f15.5,f15.5)
  endif
  write(6,1200) eps0,iter
1200 format(/5x,'eps =',1pe8.1,'  iter =',i3)
  write(6,*)
  write(6,1211) etot, etot*0.5_dp, etot*rytoev
1211 format (5x,'Etot =',f15.6,' Ry,',f15.6, ' Ha,',f15.6,' eV')
  write(6,1221) etots, etots*0.5_dp, etots*rytoev
1221 format (5x,'Etotps =',f13.6,' Ry,',f15.6,' Ha,',f15.6,' eV') 
  write(6,1231) etot-etot0
1231 format (5x,'dEtot_ae =',f15.6,' Ry') 
  write(6,1241) etots-etots0, etot-etot0 - (etots-etots0)
1241 format (5x,'dEtot_ps =',f15.6,' Ry,','   Delta E=',f15.6 )
  write(13,'(5x,''dEtot_ae ='',f15.6,'' Ry'')') etot-etot0
  write(13,'(5x,''dEtot_ps ='',f15.6,'' Ry,'',''   Delta E='', f15.6 )') &
       etots-etots0, etot-etot0-(etots-etots0)

  write(6,1251) ekin, ekin*0.5_dp, ekin*rytoev
1251 format (/,5x,'Ekin =',f15.6,' Ry,',f15.6,' Ha,',f15.6,' eV')

  write(6,1261) encl, encl*0.5_dp, encl*rytoev
1261 format (5x,'Encl =',f15.6,' Ry,',f15.6, ' Ha,',f15.6,' eV') 
  write(6,1271) ehrt, ehrt*0.5_dp, ehrt*rytoev
1271 format (5x,'Ehrt =',f15.6,' Ry,',f15.6,' Ha,',f15.6,' eV') 
  write(6,1281) ecxc, ecxc*0.5_dp, ecxc*rytoev
1281 format (5x,'Ecxc =',f15.6,' Ry,',f15.6,' Ha,',f15.6,' eV')
  if (nlcc) write(6,1282) ecc, ecc*0.5_dp, ecc*rytoev
1282 format (5x,'(Ecc =',f15.6,' Ry,',f15.6,' Ha,',f15.6,' eV)')
  write(6,1291) evxt, evxt*0.5_dp, evxt*rytoev
1291 format(5x,'Evxt =',f15.6,' Ry,',f15.6,' Ha,',f15.6,' eV')
  write(6,1292) epseu, epseu*0.5_dp, epseu*rytoev
1292 format (5x,'Epseu=',f15.6,' Ry,',f15.6, ' Ha,',f15.6,' eV') 
  if(isic.ne.0) write(6,1300) dhrsic+dxcsic, dhrsic, dxcsic
1300 format(5x,'desic:'/5x,0pf12.4,24x,2(0pf12.4))
  write(6,120)
120 format (/,5x,22('-'), ' End of pseudopotential test ',22('-'),/)
  !
  write(13,*)
  if (file_wavefunctionsps.ne.' ') then
     open(unit=16,file=file_wavefunctionsps,status='unknown', &
          err=1110, iostat=ios,form='formatted')
1110 call errore('write_resultps','opening file_wavefunctionsps',abs(ios))
     do n=1,mesh
        write(16,'(8f10.6)') r(n),(phits(n,i), &
             i=nwfts,max(1,nwfts-6),-1)
        !            write(16,'(6f12.6)') r(n),(vnl(n,i)-vpsloc(n), i=lmax,0,-1)
     enddo
     close(16)
  endif

  return
end subroutine write_resultsps
