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
  use kinds,     only : dp
  use radial_grids, only : ndmx
  use io_global, only : stdout, ionode, ionode_id
  use mp,        only : mp_bcast
  use constants, only : eps6
  use ld1inc,    only : title, rel, zed, zval, lsd, isic, latt, beta, tr2, &
                  nwfts, nnts, llts, jjts, elts, octs, iswts, enlts, nstoaets, &
                  grid, enl,  eps0, iter, etot, etots, etot0, lpaw, &
                  etots0, ekin, encl, ehrt, ecxc, nlcc, ecc, evxt, epseu, &
                  dhrsic, dxcsic, file_wavefunctionsps, phits, rytoev_fact, &
                  verbosity, frozen_core, ae_fc_energy
       
  use funct, only: get_dft_name
  implicit none

  integer :: i, j, n, m, l, ios
  real(DP) :: work(ndmx), int_0_inf_dr, ravg, sij
  character (len=20) :: dft_name
  !
  !
  dft_name = get_dft_name()
  write(stdout,110)
110 format (/,5x,22('-'),' Testing the pseudopotential ',22('-'),/)
  write(stdout,1150) title
  if(rel.eq.1) write(stdout,'(5x,''scalar relativistic calculation'')')
  if(rel.eq.2) write(stdout,'(5x,''dirac relativistic calculation'')')
1150 format(5x,a75)
  write(stdout,1250) zed, zval
1250 format(/5x,'atomic number is',f6.2,'   valence charge is',f6.2)
  write(stdout,2300) dft_name(1:len_trim(dft_name)),lsd,isic,latt,beta,tr2
2300 format(5x,'dft =',a,'   lsd =',i1,' sic =',i1,' latt =',i1, &
       '  beta=',f4.2,' tr2=',1pe7.1)
  write(stdout,1270) grid%mesh,grid%r(grid%mesh),grid%xmin,grid%dx
1270 format(5x,'mesh =',i4,' r(mesh) =',f10.5,' xmin =',f6.2,' dx =',f8.5)
  if (rel.lt.2) then
     write(stdout,1000)
1000 format(/5x,'n l     nl             e AE (Ry) ',  &
          '       e PS (Ry)    De AE-PS (Ry) ')
     do n=1,nwfts
        if (octs(n)>-eps6) write(stdout,1100) &
             nnts(n),llts(n),elts(n),iswts(n),octs(n), &
             enl(nstoaets(n)),enlts(n), &
             enl(nstoaets(n))-enlts(n)
     enddo
     if (ionode) write(13,1100)  &
          (nnts(n),llts(n),elts(n),iswts(n),octs(n), &
          enl(nstoaets(n)),enlts(n),  &
          enl(nstoaets(n))-enlts(n),  n=1,nwfts)
1100 format(4x,2i2,5x,a2,i4,'(',f5.2,')',f15.5,f15.5,f15.5)
  else
     write(stdout,1001)
1001 format(/5x,'n l  j  nl             e AE (Ry)',  &
          '       e PS (Ry)    De AE-PS (Ry) ')
     do n=1,nwfts
        if(octs(n)>-eps6) write(stdout,1101) &
             nnts(n),llts(n),jjts(n),elts(n),iswts(n),octs(n), &
             enl(nstoaets(n)),enlts(n), enl(nstoaets(n))-enlts(n)
     enddo
     if (ionode) write(13,1101)  &
          (nnts(n),llts(n),jjts(n),elts(n),iswts(n),octs(n), &
          enl(nstoaets(n)),enlts(n),  &
          enl(nstoaets(n))-enlts(n),  n=1,nwfts)
1101 format(4x,2i2,f4.1,1x,a2,i4,'(',f5.2,')',f15.5,f15.5,f15.5)
  endif
  write(stdout,1200) eps0,iter
1200 format(/5x,'eps =',1pe8.1,'  iter =',i3)
  write(stdout,*)
  write(stdout,1211) etot, etot*0.5_dp, etot*rytoev_fact
1211 format (5x,'Etot =',f15.6,' Ry,',f15.6, ' Ha,',f15.6,' eV')
  write(stdout,1221) etots, etots*0.5_dp, etots*rytoev_fact
1221 format (5x,'Etotps =',f13.6,' Ry,',f15.6,' Ha,',f15.6,' eV') 
  if (frozen_core.or.(verbosity=='high'.and.lpaw)) &
       write(stdout,1222) ae_fc_energy, ae_fc_energy*0.5_dp, &
                                           ae_fc_energy*rytoev_fact
1222 format (5x,'Etotfc =',f13.6,' Ry,',f15.6,' Ha,',f15.6,' eV') 

  if (abs(etot-etot0)> 1.d-9) then 
     write(stdout,1231) etot-etot0
1231 format (5x,'dEtot_ae =',f15.6,' Ry') 
     write(stdout,1241) etots-etots0, etot-etot0 - (etots-etots0)
1241 format (5x,'dEtot_ps =',f15.6,' Ry,','   Delta E=',f15.6,' Ry' )
     if (ionode) write(13,'(5x,''dEtot_ae ='',f15.6,'' Ry'')') etot-etot0
     if (ionode) write(13,&
       '(5x,''dEtot_ps ='',f15.6,'' Ry,'',''   Delta E='', f15.6,'' Ry'' )') &
          etots-etots0, etot-etot0-(etots-etots0)
  else
     if (ionode) write(13,1211) etot, etot*0.5_dp, etot*rytoev_fact
     if (ionode) write(13,1221) etots, etots*0.5_dp, etots*rytoev_fact
  endif
  write(stdout,1251) ekin, ekin*0.5_dp, ekin*rytoev_fact
1251 format (/,5x,'Ekin =',f15.6,' Ry,',f15.6,' Ha,',f15.6,' eV')

  write(stdout,1261) encl, encl*0.5_dp, encl*rytoev_fact
1261 format (5x,'Encl =',f15.6,' Ry,',f15.6, ' Ha,',f15.6,' eV') 
  write(stdout,1271) ehrt, ehrt*0.5_dp, ehrt*rytoev_fact
1271 format (5x,'Ehrt =',f15.6,' Ry,',f15.6,' Ha,',f15.6,' eV') 
  write(stdout,1281) ecxc, ecxc*0.5_dp, ecxc*rytoev_fact
1281 format (5x,'Ecxc =',f15.6,' Ry,',f15.6,' Ha,',f15.6,' eV')
  if (nlcc) write(stdout,1282) ecc, ecc*0.5_dp, ecc*rytoev_fact
1282 format (5x,'(Ecc =',f15.6,' Ry,',f15.6,' Ha,',f15.6,' eV)')
  if (ABS(evxt)>0.0_DP) &
     write(stdout,1291) evxt, evxt*0.5_dp, evxt*rytoev_fact
1291 format(5x,'Evxt =',f15.6,' Ry,',f15.6,' Ha,',f15.6,' eV')
  if (ABS(epseu)>0.0_DP) &
     write(stdout,1292) epseu, epseu*0.5_dp, epseu*rytoev_fact
1292 format (5x,'Epseu=',f15.6,' Ry,',f15.6, ' Ha,',f15.6,' eV') 
  if(isic.ne.0) write(stdout,1300) dhrsic+dxcsic, dhrsic, dxcsic
1300 format(5x,'desic:'/5x,0pf12.4,24x,2(0pf12.4))
  write(stdout,120)
120 format (/,5x,22('-'), ' End of pseudopotential test ',22('-'),/)
  !
  if (ionode) write(13,*)
  if (file_wavefunctionsps.ne.' ') then
     if (ionode) &
        open(unit=16,file=file_wavefunctionsps,status='unknown', &
          err=1110, iostat=ios,form='formatted')
1110 call mp_bcast(ios, ionode_id)
     call errore('write_resultps','opening file_wavefunctionsps',abs(ios))
     if (ionode) then
        do n=1,grid%mesh
           write(16,'(8f10.6)') grid%r(n),(phits(n,i), &
                                      i=nwfts,max(1,nwfts-6),-1)
        !            write(16,'(6f12.6)') r(n),(vnl(n,i)-vpsloc(n), i=lmax,0,-1)
        enddo
        close(16)
     endif
  endif

  return
end subroutine write_resultsps
