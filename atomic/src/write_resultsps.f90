!
! Copyright (C) 2004-2010 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------
SUBROUTINE write_resultsps ( )
  !--------------------------------------------------------------
  USE kinds,     ONLY : dp
  USE radial_grids, ONLY : ndmx
  USE io_global, ONLY : stdout, ionode, ionode_id
  USE mp,        ONLY : mp_bcast
  USE constants, ONLY : eps6
  USE ld1inc,    ONLY : title, rel, zed, zval, lsd, isic, latt, beta, tr2, &
                  nwfts, nnts, llts, jjts, elts, octs, iswts, enlts, nstoaets, &
                  grid, enl,  eps0, iter, etot, etots, etot0, lpaw, &
                  etots0, ekin, encl, ehrt, ecxc, nlcc, ecc, evxt, epseu, &
                  dhrsic, dxcsic, file_wavefunctionsps, phits, rytoev_fact, &
                  verbosity, frozen_core, ae_fc_energy, jj, max_out_wfc

  USE ld1inc,    ONLY : nwf, el, psi, rcut
  USE funct, ONLY: get_dft_name
  IMPLICIT NONE

  INTEGER :: counter
  real(DP) :: psiaux(ndmx,2*max_out_wfc), phase
  CHARACTER (len=2) :: elaux(2*max_out_wfc)

  INTEGER :: i, j, n, wfc_num, ios
  CHARACTER (len=20) :: dft_name
  !
  !
  dft_name = get_dft_name()
  WRITE(stdout,"(/,5x,22('-'),' Testing the pseudopotential ',22('-'),/)")
  WRITE(stdout,'(5x,a75)') title
  IF(rel==1) WRITE(stdout,'(5x,''scalar relativistic calculation'')')
  IF(rel==2) WRITE(stdout,'(5x,''dirac relativistic calculation'')')
  WRITE(stdout,"(/5x,'atomic number is',f6.2,'   valence charge is',f6.2)") &
       zed, zval
  WRITE(stdout,100) dft_name(1:len_trim(dft_name)),lsd,isic,latt,beta,tr2
100 FORMAT(5x,'dft =',a,'   lsd =',i1,' sic =',i1,' latt =',i1, &
       '  beta=',f4.2,' tr2=',1pe7.1)
  WRITE(stdout,200) grid%mesh,grid%r(grid%mesh),grid%xmin,grid%dx
200 FORMAT(5x,'mesh =',i4,' r(mesh) =',f10.5,' xmin =',f6.2,' dx =',f8.5)
  IF (rel<2) THEN
     WRITE(stdout,300)
300 FORMAT(/5x,'n l     nl             e AE (Ry) ',  &
          '       e PS (Ry)    De AE-PS (Ry) ')
     DO n=1,nwfts
        IF (verbosity=='high') THEN
           IF (octs(n)>-eps6) THEN
              IF (ABS(enl(nstoaets(n))-enlts(n))< 5.d-3) THEN
                 WRITE(stdout,401) &
                 nnts(n),llts(n),elts(n),iswts(n),octs(n), &
                 enl(nstoaets(n)),enlts(n), &
                 enl(nstoaets(n))-enlts(n)
              ELSE
!
!     put a ! close to the eigenvalues that differ more than 5 mRy
!
                 WRITE(stdout,403) &
                 nnts(n),llts(n),elts(n),iswts(n),octs(n), &
                 enl(nstoaets(n)),enlts(n), &
                 enl(nstoaets(n))-enlts(n)
              ENDIF
           ENDIF
        ELSE
           IF (octs(n)>-eps6) THEN
              IF (ABS(enl(nstoaets(n))-enlts(n))< 5.d-3) THEN
                WRITE(stdout,400) &
                  nnts(n),llts(n),elts(n),iswts(n),octs(n), &
                  enl(nstoaets(n)),enlts(n), &
                  enl(nstoaets(n))-enlts(n)
              ELSE
                 WRITE(stdout,402) &
                   nnts(n),llts(n),elts(n),iswts(n),octs(n), &
                   enl(nstoaets(n)),enlts(n), &
                   enl(nstoaets(n))-enlts(n)
              ENDIF
           ENDIF
        ENDIF
     ENDDO
     IF (ionode) WRITE(13,400)  &
          (nnts(n),llts(n),elts(n),iswts(n),octs(n), &
          enl(nstoaets(n)),enlts(n),  &
          enl(nstoaets(n))-enlts(n),  n=1,nwfts)
400 FORMAT(4x,2i2,5x,a2,i4,'(',f5.2,')',f15.5,f15.5,f15.5)
401 FORMAT(4x,2i2,5x,a2,i4,'(',f5.2,')',f15.5,f15.5,f15.8)
402 FORMAT(4x,2i2,5x,a2,i4,'(',f5.2,')',f15.5,f15.5,f15.5,"  !")
403 FORMAT(4x,2i2,5x,a2,i4,'(',f5.2,')',f15.5,f15.5,f15.8,"  !")
  ELSE
     WRITE(stdout,500)
500 FORMAT(/5x,'n l  j  nl             e AE (Ry)',  &
          '       e PS (Ry)    De AE-PS (Ry) ')
     DO n=1,nwfts
        IF (verbosity=='high') THEN
           IF(octs(n)>-eps6) THEN
             IF (ABS(enl(nstoaets(n))-enlts(n))< 5.d-3) THEN
                WRITE(stdout,601) &
                nnts(n),llts(n),jjts(n),elts(n),iswts(n),octs(n), &
                enl(nstoaets(n)),enlts(n), enl(nstoaets(n))-enlts(n)
             ELSE
                WRITE(stdout,603) &
                nnts(n),llts(n),jjts(n),elts(n),iswts(n),octs(n), &
                enl(nstoaets(n)),enlts(n), enl(nstoaets(n))-enlts(n)
             ENDIF
           ENDIF
        ELSE
           IF(octs(n)>-eps6) THEN
             IF (ABS(enl(nstoaets(n))-enlts(n))< 5.d-3) THEN
                WRITE(stdout,600) &
                  nnts(n),llts(n),jjts(n),elts(n),iswts(n),octs(n), &
                  enl(nstoaets(n)),enlts(n), enl(nstoaets(n))-enlts(n)
             ELSE
                WRITE(stdout,602) &
                   nnts(n),llts(n),jjts(n),elts(n),iswts(n),octs(n), &
                   enl(nstoaets(n)),enlts(n), enl(nstoaets(n))-enlts(n)
             ENDIF
           ENDIF
        ENDIF
     ENDDO
     IF (ionode) WRITE(13,600)  &
          (nnts(n),llts(n),jjts(n),elts(n),iswts(n),octs(n), &
          enl(nstoaets(n)),enlts(n),  &
          enl(nstoaets(n))-enlts(n),  n=1,nwfts)
600 FORMAT(4x,2i2,f4.1,1x,a2,i4,'(',f5.2,')',f15.5,f15.5,f15.5)
601 FORMAT(4x,2i2,f4.1,1x,a2,i4,'(',f5.2,')',f15.5,f15.5,f15.8)
602 FORMAT(4x,2i2,f4.1,1x,a2,i4,'(',f5.2,')',f15.5,f15.5,f15.5,"  !")
603 FORMAT(4x,2i2,f4.1,1x,a2,i4,'(',f5.2,')',f15.5,f15.5,f15.8,"  !")
  ENDIF
  WRITE(stdout,"(/5x,'eps =',1pe8.1,'  iter =',i3)") eps0,iter
  WRITE(stdout,*)
  WRITE(stdout,700) etot, etot*0.5_dp, etot*rytoev_fact
700 FORMAT (5x,'Etot =',f15.6,' Ry,',f15.6, ' Ha,',f15.6,' eV')
  WRITE(stdout,800) etots, etots*0.5_dp, etots*rytoev_fact
800 FORMAT (5x,'Etotps =',f13.6,' Ry,',f15.6,' Ha,',f15.6,' eV')
  IF (frozen_core.or.(verbosity=='high'.and.lpaw)) &
       WRITE(stdout,900) ae_fc_energy, ae_fc_energy*0.5_dp, &
                                           ae_fc_energy*rytoev_fact
900 FORMAT (5x,'Etotfc =',f13.6,' Ry,',f15.6,' Ha,',f15.6,' eV')

  IF (abs(etot-etot0)> 1.d-9) THEN
     WRITE(stdout,"(5x,'dEtot_ae =',f15.6,' Ry')") etot-etot0
     WRITE(stdout,1000) etots-etots0, etot-etot0 - (etots-etots0)
1000 FORMAT (5x,'dEtot_ps =',f15.6,' Ry,','   Delta E=',f15.6,' Ry' )
     IF (ionode) WRITE(13,'(5x,''dEtot_ae ='',f15.6,'' Ry'')') etot-etot0
     IF (ionode) WRITE(13,&
       '(5x,''dEtot_ps ='',f15.6,'' Ry,'',''   Delta E='', f15.6,'' Ry'' )') &
          etots-etots0, etot-etot0-(etots-etots0)
  ELSE
     IF (ionode) WRITE(13,700) etot, etot*0.5_dp, etot*rytoev_fact
     IF (ionode) WRITE(13,800) etots, etots*0.5_dp, etots*rytoev_fact
  ENDIF
  WRITE(stdout,1100) ekin, ekin*0.5_dp, ekin*rytoev_fact
1100 FORMAT (/,5x,'Ekin =',f15.6,' Ry,',f15.6,' Ha,',f15.6,' eV')

  WRITE(stdout,1200) encl, encl*0.5_dp, encl*rytoev_fact
1200 FORMAT (5x,'Encl =',f15.6,' Ry,',f15.6, ' Ha,',f15.6,' eV')
  WRITE(stdout,1271) ehrt, ehrt*0.5_dp, ehrt*rytoev_fact
1271 FORMAT (5x,'Ehrt =',f15.6,' Ry,',f15.6,' Ha,',f15.6,' eV')
  WRITE(stdout,1281) ecxc, ecxc*0.5_dp, ecxc*rytoev_fact
1281 FORMAT (5x,'Ecxc =',f15.6,' Ry,',f15.6,' Ha,',f15.6,' eV')
  IF (nlcc) WRITE(stdout,1282) ecc, ecc*0.5_dp, ecc*rytoev_fact
1282 FORMAT (5x,'(Ecc =',f15.6,' Ry,',f15.6,' Ha,',f15.6,' eV)')
  IF (abs(evxt)>0.0_DP) &
     WRITE(stdout,1291) evxt, evxt*0.5_dp, evxt*rytoev_fact
1291 FORMAT(5x,'Evxt =',f15.6,' Ry,',f15.6,' Ha,',f15.6,' eV')
  IF (abs(epseu)>0.0_DP) &
     WRITE(stdout,1292) epseu, epseu*0.5_dp, epseu*rytoev_fact
1292 FORMAT (5x,'Epseu=',f15.6,' Ry,',f15.6, ' Ha,',f15.6,' eV')
  IF(isic/=0) WRITE(stdout,1300) dhrsic+dxcsic, dhrsic, dxcsic
1300 FORMAT(5x,'desic:'/5x,0pf12.4,24x,2(0pf12.4))
  WRITE(stdout,120)
120 FORMAT (/,5x,22('-'), ' End of pseudopotential test ',22('-'),/)
  !
  IF (ionode) WRITE(13,*)
  !
  IF (file_wavefunctionsps/=' ') THEN
     counter=1
     wfc_num=MIN(nwfts, max_out_wfc)
     DO i=1,nwfts
        IF (counter > max_out_wfc) exit
        elaux(counter)=elts(i)
        psiaux(:,counter)=phits(:,i)
        DO j=nwf,1,-1
           IF ( elts(i) == el(j) .and. jjts(i)==jj(j) ) THEN
              DO n=grid%mesh,1,-1
                 phase = psiaux(n,counter)*psi(n,1,j)
                 IF ( abs(phase) > 1.d-12 ) THEN
                    phase = phase / abs(phase)
                    exit
                 ENDIF
              ENDDO
              psiaux(:,wfc_num+counter)=psi(:,1,j)*phase
              elaux(wfc_num+counter)=el(j)
              exit
           ENDIF
        ENDDO
        counter=counter+1
     ENDDO
     counter = counter - 1
     CALL write_wfcfile(file_wavefunctionsps,psiaux,elaux,2*counter)
  ENDIF

  RETURN
END SUBROUTINE write_resultsps
