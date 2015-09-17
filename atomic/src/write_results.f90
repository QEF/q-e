!
! Copyright (C) 2004-2012 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------
subroutine write_results 
  !--------------------------------------------------------------
  use radial_grids, only : ndmx
  use kinds,        only : dp
  use io_global, only : stdout, ionode_id, ionode
  use mp,        only : mp_bcast
  use mp_world,  only : world_comm
  use constants, only : eps6
  use ld1inc,    only : title, rel, zed, lsd, nspin, isic, latt, beta, tr2, &
                        grid, enzero, etot, ekin, encl, ehrt, evxt, ecxc, &
                        ehrtcc, ehrtcv, ehrtvv, enclv, enclc, verbosity,  &
                        nwf, nn, ll, jj, el, isw, oc, enl, file_wavefunctions, &
                        dhrsic, dxcsic, eps0, iter, psi, rytoev_fact, lsmall, &
                        core_state, ekinc, ekinv, ae_fc_energy, cau_fact, &
                        relpert, evel, edar, eso, noscf, iswitch, rho, &
                        file_charge, max_out_wfc

  use funct, only :  get_iexch, get_dft_name, write_dft_name
  implicit none

  integer :: is, i, j, n, m, im(40), ios, counter, ismax
  real(DP):: work(ndmx), dum, int_0_inf_dr, ravg, r2avg, sij, ene, mm, &
             sij1, sij2, charge_large, charge_small, work1(ndmx), work2(ndmx)
  real(DP) :: psiaux(ndmx,max_out_wfc)
  logical :: ok, oep, print_fc
  character (len=20) :: dft_name
  character (len=2) :: elaux(max_out_wfc)
  character (len=60) :: vstates
  character (len=256) :: nomefile
  character (len=6), dimension(2) :: suffix
  !
  !
  dft_name = get_dft_name()
  write(stdout,"(5x,27('-'),' All-electron run ',28('-'),/)")
  write(stdout,1150) title
  if(rel.eq.1) write(stdout,"(5x,'scalar relativistic calculation')")
  if(rel.eq.2) write(stdout,"(5x,'dirac relativistic calculation')")
1150 format(5x,a75)
  if (zed.ne.0.0) write(stdout,1250) zed
1250 format(/5x,'atomic number is',f6.2)
  write(stdout,2300) dft_name(1:len_trim(dft_name)),lsd,isic,latt,beta,tr2
  CALL write_dft_name () 
2300 format(5x,'dft =',a,'   lsd =',i1,' sic =',i1,' latt =',i1, &
       '  beta=',f4.2,' tr2=',1pe7.1)
  write(stdout,1270) grid%mesh,grid%r(grid%mesh),grid%xmin,grid%dx
1270 format(5x,'mesh =',i4,' r(mesh) =',f10.5,' a.u. xmin =',f6.2,' dx =',f8.5)
  if (rel==0 .and. .not.relpert) then
     write(stdout,'(5x,"1 Ry = ",f12.8, " eV" )') rytoev_fact
  else
     write(stdout,'(5x,"1 Ry = ",f12.8, " eV, c = ",f12.8 )') rytoev_fact, &
                                                            cau_fact
  endif
  if (rel.lt.2) then
     write(stdout,1000)
1000 format(/5x, &
          'n l     nl                  e(Ry) ','         e(Ha)          e(eV)')

     oep = get_iexch() .eq. 4
     if (oep) enl(1:nwf) = enl(1:nwf) - enzero(isw(1:nwf))
     do n=1,nwf
        if (oc(n)>-eps6) then 
          IF (verbosity=='high') THEN
             write(stdout,1103) &
                nn(n),ll(n),el(n),isw(n),oc(n),enl(n),enl(n)*0.5_dp, &
                enl(n)*rytoev_fact
          ELSE
             write(stdout,1100) &
                nn(n),ll(n),el(n),isw(n),oc(n),enl(n),enl(n)*0.5_dp, &
                enl(n)*rytoev_fact
          ENDIF
          !!! relativistic perturbative terms
          if (verbosity=='high'.and.relpert) then
             if(rel.eq.0) write(stdout,1102) &
                "-E_vel",evel(n),evel(n)*0.5_dp,evel(n)*rytoev_fact
             if(rel.eq.0) write(stdout,1102) &
                "-E_Dar",edar(n),edar(n)*0.5_dp,edar(n)*rytoev_fact
             if(rel.eq.0 .and. ll(n).gt.0) write(stdout,1102) &
                "-E_S-O",eso(n),eso(n)*0.5_dp,eso(n)*rytoev_fact
          endif
          !!!
        endif     
     enddo
     !!!
1102 format(18x,a6,f15.4,f15.4,f15.4)
     !!!

     if (oep) then
        enl(1:nwf) = enl(1:nwf) + enzero(isw(1:nwf))
        write(stdout,*) 
!!1100 format(4x,2i2,5x,a2,i2,'(',f5.2,')',f15.4,f15.4,f15.4)
        write(stdout,'(5x,a)') "OEP WARNING: printed eigenvalues were shifted by"
        if (nspin==1) write(stdout,'(17x,a,3f15.4)') ( "shift :", &
                            enzero(is), enzero(is)*0.5_dp, &
                            enzero(is)*rytoev_fact, is=1,nspin)
        if (nspin==2) write(stdout,'(8x,a,i2,3x,a,3f15.4)') ( "spin",is,"shift :", &
                            enzero(is), enzero(is)*0.5_dp, &
                            enzero(is)*rytoev_fact, is=1,nspin)
     end if
  else
     write(stdout,1001)
1001 format(/5x, &
          'n l j   nl                  e(Ry) ','         e(Ha)          e(eV)')
     write(stdout,'(5x,"Spin orbit split results")')
     do n=1,nwf
        IF (verbosity=='high') THEN
           if (oc(n)>-eps6) write(stdout,1123) &
               nn(n),ll(n),jj(n),el(n),isw(n),oc(n),enl(n),enl(n)*0.5_dp, &
               enl(n)*rytoev_fact
        ELSE
           if (oc(n)>-eps6) write(stdout,1120) &
              nn(n),ll(n),jj(n),el(n),isw(n),oc(n),enl(n),enl(n)*0.5_dp, &
              enl(n)*rytoev_fact
        ENDIF
     enddo
     write(stdout,'(5x,"Averaged results")')
     ok=.true.
     do n=1,nwf
        if (ll(n).gt.0.and.ok) then
           if (oc(n)+oc(n+1)>-eps6) then
              ene=(enl(n)*2.0_dp*ll(n) &
                   + enl(n+1)*(2.0_dp*ll(n)+2.0_dp))/(4.0_dp*ll(n)+2.0_dp)
              IF (verbosity=='high') THEN
                 write(stdout,1103) nn(n),ll(n),el(n),isw(n),oc(n)+oc(n+1), &
                      ene, ene*0.5_dp, ene*rytoev_fact
              ELSE
                 write(stdout,1100) nn(n),ll(n),el(n),isw(n),oc(n)+oc(n+1), &
                      ene, ene*0.5_dp, ene*rytoev_fact
              ENDIF
              ok=.false.
           endif
        elseif (ll(n)==0) then
           if ( oc(n) > -eps6 ) then
              IF (verbosity=='high') THEN
                 write(stdout,1103) nn(n),ll(n),el(n),isw(n),oc(n), &
                       enl(n), enl(n)*0.5_dp, enl(n)*rytoev_fact
              ELSE
                 write(stdout,1100) nn(n),ll(n),el(n),isw(n),oc(n), &
                       enl(n), enl(n)*0.5_dp, enl(n)*rytoev_fact
              ENDIF
              ok=.true.
           endif
        else
           ok=.true.
        endif
     enddo
  endif

  !!!
  !!! eigenvalues with perturbative relativistic corrections
  !!!
  if ( relpert ) then
    !
     write(stdout,'(/,5x,a)') &
        "relativistic eigenvalues in first order perturbation theory"
     write(stdout,1001)
     do n=1,nwf
       if (oc(n)>-eps6) then
         if (ll(n).gt.0) then
           ! j = l-1/2
           ene = enl(n) - evel(n) - edar(n) + eso(n)*real(ll(n)+1,DP)
           IF (verbosity=='high') THEN
              write(stdout,1124) &
                nn(n),ll(n),ll(n)-0.5_DP,el(n),isw(n),&
                ene,ene*0.5_dp, ene*rytoev_fact
           ELSE
              write(stdout,1121) &
                nn(n),ll(n),ll(n)-0.5_DP,el(n),isw(n),&
                ene,ene*0.5_dp, ene*rytoev_fact
           ENDIF
           ! j = l+1/2
           ene = enl(n) - evel(n) - edar(n) - eso(n)*real(ll(n),DP)
           IF (verbosity=='high') THEN
              write(stdout,1124) &
                 nn(n),ll(n),ll(n)+0.5_DP,el(n),isw(n),&
                 ene,ene*0.5_dp, ene*rytoev_fact
           ELSE
              write(stdout,1121) &
                 nn(n),ll(n),ll(n)+0.5_DP,el(n),isw(n),&
                 ene,ene*0.5_dp, ene*rytoev_fact
           ENDIF
            ! avrebbe senso lasciare le occupazioni?
         else
           ! j = l+1/2
           ene = enl(n) - evel(n) - edar(n)
           IF (verbosity=='high') THEN
              write(stdout,1124) &
                 nn(n),ll(n),ll(n)+0.5_DP,el(n),isw(n),&
                 ene,ene*0.5_dp, ene*rytoev_fact
           ELSE
               write(stdout,1121) &
                 nn(n),ll(n),ll(n)+0.5_DP,el(n),isw(n),&
                 ene,ene*0.5_dp, ene*rytoev_fact
           ENDIF
         endif
       endif
     enddo
    !
  endif
1121 format(4x,2i2,f4.1,1x,a2,i2,7x,f15.4,f15.4,f15.4)
1124 format(4x,2i2,f4.1,1x,a2,i2,7x,f18.9,f18.9,f15.7)
  !!!
1100 format(4x,2i2,5x,a2,i2,'(',f5.2,')',f15.4,f15.4,f15.4)
1120 format(4x,2i2,f4.1,1x,a2,i2,'(',f5.2,')',f15.4,f15.4,f15.4)
1103 format(4x,2i2,5x,a2,i2,'(',f5.2,')',f18.9,f18.9,f15.7)
1123 format(4x,2i2,f4.1,1x,a2,i2,'(',f5.2,')',f18.9,f18.9,f15.7)


  if (noscf) goto 500

  write(stdout,1200) eps0, iter
1200 format(/5x,'final scf error: ',1pe8.1,' reached in ',i3,' iterations')
  print_fc=(verbosity=='high'.and.iswitch>1) 
  write(stdout,*)
  if (print_fc) then
      vstates=' '
      do n=1,nwf
         if (.not.core_state(n)) vstates=TRIM(vstates)//el(n)//","
      enddo
      write(stdout,'(5x,"The valence states are: ", a60,/)') vstates 
  endif    
  if (print_fc) write(6,"(5x,'Total energy of the atom:',/)")
  write(stdout,"(5x,'Etot =',f15.6,' Ry,',f15.6,' Ha,',f15.6,' eV')") &
       etot, etot*0.5_dp, etot*rytoev_fact
  if (lsd.eq.1) then
     mm=0.d0
     do n=1,nwf
        if (oc(n).gt.0.0_dp) then
           if (isw(n).eq.1) mm=mm+oc(n)
           if (isw(n).eq.2) mm=mm-oc(n)
        endif
     enddo
     write(stdout,"(5x,'Total magnetization:',f8.2,' Bohr mag. ')") mm
  endif
  if (print_fc) write(stdout,'(/,5x,"Kinetic energy:")')
  write(stdout,"(/,5x,'Ekin =',f15.6,' Ry,',f15.6,' Ha,',f15.6,' eV')")&
       ekin, ekin*0.5_dp,  ekin*rytoev_fact
  if (print_fc) then
     write(stdout,"(5x,'Ekinc=',f15.6,' Ry,',f15.6,' Ha,',f15.6,' eV')")&
                   ekinc, ekinc*0.5_dp,  ekinc*rytoev_fact
     write(stdout,"(5x,'Ekinv=',f15.6,' Ry,',f15.6,' Ha,',f15.6,' eV')")&
                   ekinv, ekinv*0.5_dp,  ekinv*rytoev_fact
     write(stdout,*)
  endif
  if (print_fc) write(6,'(/,5x,"Interaction between the nucleus and &
                                       &the electrons:",/)')
  write(stdout,"(5x,'Encl =',f15.6,' Ry,',f15.6,' Ha,',f15.6,' eV')")&
       encl, encl*0.5_dp, encl*rytoev_fact
  if (print_fc) then
     write(stdout, "(5x,'Enclc=',f15.6,' Ry,',f15.6,' Ha,',f15.6,' eV')") &
       enclc, enclc*0.5_dp, enclc*rytoev_fact
     write(stdout, "(5x,'Enclv=',f15.6,' Ry,',f15.6,' Ha,',f15.6,' eV')") &
       enclv, enclv*0.5_dp, enclv*rytoev_fact
     write(stdout,*)
  endif   
  if (print_fc) write(6,'(/,5x,"Hartree energy:",/)')
  write(stdout,"(5x,'Eh   =',f15.6,' Ry,',f15.6, ' Ha,',f15.6,' eV')") &
       ehrt, ehrt*0.5_dp, ehrt*rytoev_fact
  if (print_fc) then
     write(stdout, "(5x,'Ehcc =',f15.6,' Ry,',f15.6,' Ha,',f15.6,' eV')") &
       ehrtcc, ehrtcc*0.5_dp, ehrtcc*rytoev_fact
     write(stdout, "(5x,'Ehcv =',f15.6,' Ry,',f15.6,' Ha,',f15.6,' eV')") &
       ehrtcv, ehrtcv*0.5_dp, ehrtcv*rytoev_fact
     write(stdout,"(5x,'Ehvv =',f15.6,' Ry,',f15.6,' Ha,',f15.6,' eV')") &
       ehrtvv, ehrtvv*0.5_dp, ehrtvv*rytoev_fact
     write(stdout,*)
  endif
  if (print_fc) write(stdout, "(/,5x,'Exchange and correlation energy:',/)") 
  write(stdout,"(5x,'Exc  =',f15.6,' Ry,',f15.6,' Ha,',f15.6,' eV')") &
       ecxc, ecxc*0.5_dp, ecxc*rytoev_fact
  if (print_fc)  write(stdout,*)
  if (ABS(evxt)>0.0_DP) then
     if (verbosity=='high') &
          write(stdout,"(/,5x,'Interaction with the external potential:')") 
     write(stdout, "(5x,'Evxt =',f15.6,' Ry,',f15.6,' Ha,',f15.6,' eV')") &
       evxt, evxt*0.5_dp, evxt*rytoev_fact
  endif
  if (print_fc) then
     write(stdout, "(/,5x,'Estimated frozen-core energy from all-electron calculation:')") 
     write(stdout, "(5x,'Efc = Ekinv + Enclv + Ehvv + Ehcv + Exc')") 
     write(stdout, "(5x,'Ed = Etot - Efc',/)") 
     write(stdout, "(5x,'Efc  =',f15.6,' Ry,',f15.6,' Ha,',f15.6,' eV')") &
       ae_fc_energy, ae_fc_energy*0.5_dp, ae_fc_energy*rytoev_fact
     write(stdout, "(5x,'Ed   =',f15.6,' Ry,',f15.6,' Ha,',f15.6,' eV')") &
      etot-ae_fc_energy, (etot-ae_fc_energy)*0.5_dp, (etot-ae_fc_energy)*rytoev_fact
  endif

  if (isic.ne.0) then
     write(stdout,*)
     write(stdout,'(5x,"SIC information:")') 
     write(stdout,1300) dhrsic, dhrsic*0.5_dp, dhrsic*rytoev_fact
     write(stdout,2310) dxcsic, dxcsic*0.5_dp, dxcsic*rytoev_fact
     write(stdout,2320) dxcsic+dhrsic,(dxcsic+dhrsic)*0.5_dp, &
                       (dxcsic+dhrsic)*rytoev_fact  
     write(stdout,*)
     write(stdout,2311) ecxc-dxcsic-dhrsic, &
          &  (ecxc-dxcsic-dhrsic)*0.5_dp, (ecxc-dxcsic-dhrsic)*rytoev_fact 
     write(stdout,2312) ecxc-dhrsic, &
          &               (ecxc-dhrsic)*0.5_dp, (ecxc-dhrsic)*rytoev_fact 
     write(stdout,2313) ehrt+dhrsic, &
          &              (ehrt+dhrsic)*0.5_dp, (ehrt+dhrsic)*rytoev_fact
1300 format(5x,'Esich=',f15.6,' Ry,',f15.6,' Ha,',f15.6,' eV') 
2310 format(5x,'Esicxc=',f14.6,' Ry,',f15.6,' Ha,',f15.6,' eV') 
2311 format(5x,'tot-Exc=',f13.6,' Ry,',f15.6,' Ha,',f15.6,' eV') 
2312 format(5x,'int-Exc=',f13.6,' Ry,',f15.6,' Ha,',f15.6,' eV') 
2313 format(5x,'int-Eh=',f14.6,' Ry,',f15.6,' Ha,',f15.6,' eV') 
2320 format(5x,'Esictot=',f13.6,' Ry,',f15.6,' Ha,',f15.6,' eV') 
  endif

500 continue

  write(stdout,1310)
1310 format(//5x,'normalization and overlap integrals'/)

  do i=1,nwf
     dum=0.0_dp
     do m=1,grid%mesh
        dum=max(dum,abs(psi(m,1,i)))
        if(dum.eq.abs(psi(m,1,i)))im(i)=m
     enddo
  enddo

  charge_large=0.0_DP
  charge_small=0.0_DP
  do i=1,nwf
     do j=i,nwf
        if (ll(i)==ll(j).and.jj(i)==jj(j).and.isw(i).eq.isw(j).and. &
            oc(i).ge.-1.d-12.and.oc(j).ge.-1.d-12) then
           if (rel<2) then
              do m=1,grid%mesh
                 work(m)=psi(m,1,i)*psi(m,1,j)
              enddo
           else
              do m=1,grid%mesh
                 work(m)=psi(m,1,i)*psi(m,1,j)+psi(m,2,i)*psi(m,2,j)
              enddo
              if (i==j) then
                 do m=1,grid%mesh
                    work1(m)=psi(m,1,i)*psi(m,1,j)
                    work2(m)=psi(m,2,i)*psi(m,2,j)
                 enddo
              endif
           endif
           sij = int_0_inf_dr(work,grid,grid%mesh,2*ll(i)+2)
           if (i.eq.j) then
              do m=1,grid%mesh
                 work(m)=work(m)*grid%r(m)
              enddo
              ravg = int_0_inf_dr(work,grid,grid%mesh,2*ll(i)+3)
              do m=1,grid%mesh
                 work(m)=work(m)*grid%r(m)
              enddo
              r2avg = int_0_inf_dr(work,grid,grid%mesh,2*ll(i)+4)
              write(stdout,1400) el(i),el(j),sij, ravg, r2avg, grid%r(im(i))
              IF (rel==2.and.verbosity=='high') THEN
                 sij1 = int_0_inf_dr(work1,grid,grid%mesh,2*ll(i)+2)
                 sij2 = int_0_inf_dr(work2,grid,grid%mesh,2*ll(i)+2)
                 WRITE(stdout,'(5x,"LC norm =",f12.8,&
                            &" + SC norm =",f12.8," =",f12.8 )') sij1, sij2,&
                               sij1+sij2
                 charge_large=charge_large + oc(i)*sij1
                 charge_small=charge_small + oc(i)*sij2
              ENDIF
           else
              write(stdout,1401) el(i),el(j),sij
           endif
        endif
     enddo
  enddo
1400 format(5x,'s(',a2,'/',a2,') =',f10.6,2x, &
       '<r> =',f9.4,2x,'<r2> =',f10.4,2x,'r(max) =',f9.4)
1401 format(5x,'s(',a2,'/',a2,') =',f10.6)

  IF (rel==2.and.verbosity=='high'.and..not.noscf) &
     write(stdout,'(/,5x,"LC charge =",f12.8," + SC charge =",&
            & f12.8," =",f12.8)') charge_large, charge_small, &
                                  charge_large+charge_small

  if (file_wavefunctions.ne.' ') then
     nomefile=TRIM(file_wavefunctions)
     suffix(1) = ' '
     ismax=1
     if (rel == 2 .and. lsmall ) then
        suffix(2) = '.small'
        ismax=2
     else if ( rel < 2 .and. nspin == 2 ) then
        suffix(1) = '.up'
        suffix(2) = '.dw'
        ismax=2
     end if
     do is=1,ismax
        nomefile=TRIM(file_wavefunctions)//TRIM(suffix(is))
        counter=1
        if (  rel < 2 )  then
           do i=nwf,1,-1
              if ( isw(i)==is ) then
                 elaux(counter)=el(i)
                 psiaux(:,counter)=psi(:,1,i)
                 counter=counter+1
                 if (counter>max_out_wfc) exit
              end if
           enddo
        else 
           do i=nwf,1,-1
              elaux(counter)=el(i)
              psiaux(:,counter)=psi(:,is,i)
              counter=counter+1
              if (counter>max_out_wfc) exit
           enddo
        endif
        call write_wfcfile(nomefile,psiaux,elaux,counter-1)
     enddo
  endif

  if (file_charge.ne.' ') then
     if (ionode) &
        open(unit=20,file=file_charge, status='unknown', iostat=ios, err=100 )
100  call mp_bcast(ios, ionode_id, world_comm)
     call errore('write_results','opening file'//file_charge,abs(ios))
     if (ionode) then
        IF (lsd<1) THEN
           do n=1,grid%mesh
              write(20,'(2e20.9)') grid%r(n), rho(n,1)
           enddo
        ELSE
           do n=1,grid%mesh
              write(20,'(3e20.9)') grid%r(n), rho(n,1), rho(n,2)
           enddo
        ENDIF
        close(20)
     endif
  endif

  write(stdout,"(/,5x,24('-'), ' End of All-electron run ',24('-'),/)")

  return
end subroutine write_results
