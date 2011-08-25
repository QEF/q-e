!
! Copyright (C) 2001-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE ep_matrix_element_wannier()
  !-----------------------------------------------------------------------
  !
  ! Electron-phonon calculation from data saved in fildvscf
  !
  USE kinds, ONLY : DP
  USE cell_base, ONLY : celldm, omega, ibrav
  USE ions_base, ONLY : nat, ntyp => nsp, ityp, tau, pmass
  USE gvecs, ONLY: doublegrid
  USE fft_base, ONLY : dfftp, dffts
  USE noncollin_module, ONLY : nspin_mag
  USE dynmat, ONLY : dyn, w2
  USE qpoint, ONLY : xq
  USE modes,  ONLY : npert, nirr
  USE control_ph, ONLY : trans, elph_mat
  USE units_ph, ONLY : iudyn, lrdrho, iudvscf
  USE io_global, ONLY : stdout
  USE modes, ONLY : u
  !
  IMPLICIT NONE
  !
  LOGICAL :: read_dvscf_cart, ascii_dvscf
  INTEGER :: irr, imode0, ipert, is
  ! counter on the representations
  ! counter on the modes
  ! the change of Vscf due to perturbations
  COMPLEX(DP), POINTER :: dvscfin(:,:,:), dvscfins (:,:,:)


  CALL start_clock ('elphon')

  ascii_dvscf=.false.
  if(elph_mat) read_dvscf_cart=.true.

   if(read_dvscf_cart) then
      write(stdout,*)
      write(stdout,*) 'Reading dvscf in cartesian coordinates !'
      write(stdout,*)

      u=CMPLX(0.d0,0.d0)
      do irr=1,3*nat
         u(irr,irr)=CMPLX(1.d0,0.d0)
      enddo

!      if(ascii_dvscf) then
!         ALLOCATE (dvrot ( nrxx , nspin , 3*nat) )
!         fildvscf_asc=trim(tmp_dir)//trim(prefix)//"."//trim(fildvscf)//'1'
!         open(unit=7899,file=fildvscf_asc,status='unknown')
!         dvrot=CMPLX(0.0,0.0)
!         do na=1,nat
!            do ipol=1,3
!               irr=(na-1)*3+ipol
!               do  k = 1, dfftp%nr3
!                  do j = 1, dfftp%nr2
!                     do i = 1, dfftp%nr1
!                        read(7899,*)   n, rep,imp
!                        dvrot(n,1,irr)=CMPLX(rep,imp)
!                     enddo
!                  enddo
!               enddo
!            enddo
!         enddo
!         close(7899)
!      endif

   endif

   

  !
  ! read Delta Vscf and calculate electron-phonon coefficients
  !
  imode0 = 0
  DO irr = 1, nirr
     ALLOCATE (dvscfin (dfftp%nnr, nspin_mag , npert(irr)) )
!     if(ascii_dvscf) then
!        DO ipert = 1, npert(irr)
!           dvscfin(1:dfftp%nnr,:,ipert)=dvrot(1:dfftp%nnr,:,imode0+ipert)
!        enddo
!     else
     DO ipert = 1, npert (irr)
        CALL davcio_drho ( dvscfin(1,1,ipert),  lrdrho, iudvscf, &
             imode0 + ipert,  -1 )
     END DO
!     endif
     IF (doublegrid) THEN
        ALLOCATE (dvscfins (dffts%nnr, nspin_mag , npert(irr)) )
        DO is = 1, nspin_mag
           DO ipert = 1, npert(irr)
              CALL cinterpolate (dvscfin(1,is,ipert),dvscfins(1,is,ipert),-1)
           ENDDO
        ENDDO
     ELSE
        dvscfins => dvscfin
     ENDIF
     CALL newdq (dvscfin, npert(irr))
     CALL elphel (npert (irr), imode0, dvscfins)
     !
     imode0 = imode0 + npert (irr)
     IF (doublegrid) DEALLOCATE (dvscfins)
     DEALLOCATE (dvscfin)
  ENDDO
  !
  ! now read the eigenvalues and eigenvectors of the dynamical matrix
  ! calculated in a previous run
  !
  IF (.NOT.trans) CALL readmat (iudyn, ibrav, celldm, nat, ntyp, &
       ityp, omega, pmass, tau, xq, w2, dyn)
  !
  CALL stop_clock ('elphon')
  RETURN
END SUBROUTINE ep_matrix_element_wannier

!-----------------------------------------------------------------------
SUBROUTINE elphsum_wannier
  !-----------------------------------------------------------------------
  !
  !      Sum over BZ of the electron-phonon matrix elements el_ph_mat
  !      Original routine written by Francesco Mauri
  !      Adapted to wannier functions by Matteo Calandra
  !            missing calc_sigma_yet
  !-----------------------------------------------------------------------
  USE kinds, ONLY : DP
  USE constants, ONLY : pi
  USE ions_base, ONLY : nat, ityp, tau,amass,tau, ntyp => nsp, atm
  USE cell_base, ONLY : at, bg, ibrav, celldm 
  USE fft_base,  ONLY: dfftp
  USE symm_base, ONLY : s, sr, irt, nsym, time_reversal, invs
  USE klist, ONLY : xk, nelec, nks, wk
  USE wvfct, ONLY : nbnd, et
  USE el_phon
  USE mp_global, ONLY : me_pool, root_pool, inter_pool_comm, npool, intra_pool_comm
  USE io_global, ONLY : stdout
  USE klist, only : degauss,ngauss
  USE control_flags, ONLY : modenum, noinv
  USE units_ph,       ONLY :iudyn
  USE io_files,  ONLY : prefix
  USE qpoint, ONLY : xq, nksq
  USE dynmat, ONLY : dyn, w2
  USE modes, ONLY : u,rtau, irgq, nsymq,irotmq, minus_q
  USE control_ph, only : lgamma
  USE lsda_mod, only : isk,nspin, current_spin,lsda
  USE mp,        ONLY: mp_sum
  !
  IMPLICIT NONE
  ! eps = 20 cm^-1, in Ry
  REAL(DP) :: eps
  PARAMETER (eps = 20.d0 / 13.6058d0 / 8065.5d0)
  !
  !
  INTEGER :: ik, ikk, ikq, isig, ibnd, jbnd, ipert, jpert, nu, mu, &
       vu, ngauss1, nsig, iuelph, ios, iuelphmat,icnt,i,j,rrho,nt,k
  INTEGER :: na,nb,icar,jcar,iu_sym,nmodes
  INTEGER :: iu_Delta_dyn,iu_analdyn,iu_nonanaldyn
  INTEGER :: io_file_unit
  !   for star_q
  INTEGER :: nsymloc, sloc(3,3,48), invsloc(48), irtloc(48,nat), &
             nqloc, isqloc(48), imqloc
  REAL(DP) :: rtauloc(3,48,nat), sxqloc(3,48)
  !   end of star_q definitions
  REAL(DP) :: weight, w0g1, w0g2, w0gauss, wgauss,degauss1, dosef, &
       ef1, phase_space, lambda, gamma, wg1, w0g,wgp,deltae
  REAL(DP), EXTERNAL :: dos_ef, efermig
  REAL(DP) xk_dummy(3)
  COMPLEX(DP), allocatable :: phi(:,:,:,:),phi_nonanal(:,:,:,:)
  COMPLEX(DP), allocatable :: dyn_mat_r(:,:),zz(:,:)
  CHARACTER(len=20) :: char_deg
  CHARACTER(len=1) :: char_ng
  character(len=80) :: filelph
  CHARACTER(len=20) ::  file_elphmat
  !
  COMPLEX(DP) :: el_ph_sum (3*nat,3*nat), dyn_corr(3*nat,3*nat)

  INTEGER, EXTERNAL :: find_free_unit

  nmodes=3*nat

  write(filelph,'(A5,f9.6,A1,f9.6,A1,f9.6)') 'elph.',xq(1),'.',xq(2),'.',xq(3)
  write(file_elphmat,'(A8)') 'elph.mat'
  ! parallel case: only first node writes
  IF ( me_pool /= root_pool ) THEN
     iuelph = 0
  ELSE
     !
     iuelph = find_free_unit()
     OPEN (unit = iuelph, file = filelph, status = 'unknown', err = &
          100, iostat = ios)
100  CALL errore ('elphon', 'opening file '//filelph, ABS (ios) )
     REWIND (iuelph)
     !
  END IF

  WRITE (iuelph, '(3f15.8,2i8)') xq, nsig, 3 * nat
  WRITE (iuelph, '(6e14.6)') (w2 (nu) , nu = 1, nmodes)
  
  ! parallel case: only first node writes
  IF ( me_pool /= root_pool ) THEN
     iuelphmat = 0
  ELSE
     !
     iuelphmat = find_free_unit()
     OPEN (unit = iuelphmat, file = file_elphmat, status = 'unknown', err = &
          111, iostat = ios, form='unformatted')
111  CALL errore ('elphon', 'opening file'//file_elphmat, ABS (ios) )
     REWIND (iuelphmat)
     xk_dummy(:)=xq(:)
     call cryst_to_cart(1,xk_dummy,at,-1)
     WRITE (iuelphmat) xk_dummy
     WRITE (iuelphmat) nelec
     WRITE (iuelphmat) elph_nbnd_min,elph_nbnd_max,nbnd
     WRITE (iuelphmat) nmodes, nksq, nat, ntyp
     WRITE (iuelphmat) ibrav,(celldm(j),j=1,6)
     WRITE(iuelphmat)  (atm(j),j=1,ntyp),(amass(j),j=1,ntyp), &
             (ityp(j),j=1,nat),((tau(j,i),j=1,3),i=1,nat)
     WRITE (iuelphmat) (w2 (nu) , nu = 1, nmodes)
     WRITE (iuelphmat) ((u(ipert,jpert),ipert=1,nmodes),jpert=1,nmodes)
     WRITE (iuelphmat) ((dyn(ipert,jpert),ipert=1,3*nat),jpert=1,3*nat)
     do ik=1,nksq
        IF (lgamma) THEN
           ikk = ik
           ikq = ik
        ELSE
           ikk = 2 * ik - 1
           ikq = ikk + 1
        ENDIF
        xk_dummy(:)=xk(:,ikk)
        call cryst_to_cart(1,xk_dummy,at,-1)
        WRITE (iuelphmat) (xk_dummy(ipert),ipert=1,3)
        WRITE (iuelphmat) (et(ibnd,ikk),ibnd=elph_nbnd_min,elph_nbnd_max)
        do nu=1,nmodes
           WRITE (iuelphmat) &
                ((el_ph_mat (jbnd, ibnd, ik, nu),jbnd=elph_nbnd_min,elph_nbnd_max),&
                ibnd=elph_nbnd_min,elph_nbnd_max)
        enddo
     enddo
     
     
     close(iuelphmat)
     
     !
  END IF
  !

  DO isig = 1, el_ph_nsigma
     !     degauss1 = 0.01 * isig
     degauss1 = el_ph_sigma * isig
     write(stdout,*) degauss1
     el_ph_sum(:,:) = (0.d0, 0.d0)
     phase_space = 0.d0
     !
     ! Recalculate the Fermi energy Ef=ef1 and the DOS at Ef, dosef = N(Ef)
     ! for this gaussian broadening
     !
     ! Note that the weights of k+q points must be set to zero for the
     ! following call to yield correct results
     !
      
     ef1 = efermig (et, nbnd, nks, nelec, wk, degauss1, el_ph_ngauss, 0, isk)
     dosef = dos_ef (el_ph_ngauss, degauss1, ef1, et, wk, nks, nbnd)
     ! N(Ef) is the DOS per spin, not summed over spin
     dosef = dosef / 2.d0
     !
     ! Sum over bands with gaussian weights
     !
     
     DO ik = 1, nksq
        
        !
        ! see subroutine elphel for the logic of indices
        !
        IF (lgamma) THEN
           ikk = ik
           ikq = ik
        ELSE
           ikk = 2 * ik - 1
           ikq = ikk + 1
        ENDIF
        DO ibnd = 1, nbnd
           w0g1 = w0gauss ( (ef1 - et (ibnd, ikk) ) / degauss1, ngauss1) &
                / degauss1
           xk_dummy(:)=xk(:,ikk)
           call cryst_to_cart(1,xk_dummy,at,-1)
           DO jbnd = 1, nbnd
              w0g2 = w0gauss ( (ef1 - et (jbnd, ikq) ) / degauss1, ngauss1) &
                   / degauss1
              ! note that wk(ikq)=wk(ikk)
              weight = wk (ikk) * w0g1 * w0g2
              DO jpert = 1, 3 * nat
                 DO ipert = 1, 3 * nat
                    el_ph_sum (ipert, jpert) = el_ph_sum (ipert, jpert)  +  weight * &
                         CONJG (el_ph_mat (jbnd, ibnd, ik, ipert) ) * &
                         el_ph_mat (jbnd, ibnd, ik, jpert)
                 ENDDO
              ENDDO
              phase_space = phase_space+weight
           ENDDO
        ENDDO
        
     ENDDO
     
     ! el_ph_sum(mu,nu)=\sum_k\sum_{i,j}[ <psi_{k+q,j}|dvscf_q(mu)*psi_{k,i}>
     !                                  x <psi_{k+q,j}|dvscf_q(nu)*psi_{k,i}>
     !                                  x \delta(e_{k,i}-Ef) \delta(e_{k+q,j}
     !
     ! collect contributions from all pools (sum over k-points)
     !
     
!     CALL poolreduce (2 * 3 * nat * 3 * nat, el_ph_sum)
!     CALL poolreduce (1, phase_space)
     call mp_sum ( el_ph_sum , inter_pool_comm )
     call mp_sum ( phase_space , inter_pool_comm )

     !
     ! symmetrize el_ph_sum(mu,nu) : it transforms as the dynamical matrix
     !
     
     CALL symdyn_munu (el_ph_sum, u, xq, s, invs, rtau, irt, irgq, at, &
          bg, nsymq, nat, irotmq, minus_q)
     !
     WRITE (6, 9000) degauss1, ngauss1
     WRITE (6, 9005) dosef, ef1 * 13.6058
     WRITE (6, 9006) phase_space
     IF (iuelph.NE.0) THEN
        WRITE (iuelph, 9000) degauss1, ngauss1
        WRITE (iuelph, 9005) dosef, ef1 * 13.6058
     ENDIF
     
     DO nu = 1, nmodes
        gamma = 0.0
        DO mu = 1, 3 * nat
           DO vu = 1, 3 * nat
              gamma = gamma + DBLE (CONJG (dyn (mu, nu) ) * el_ph_sum (mu, vu)&
                   * dyn (vu, nu) )
           ENDDO
        ENDDO
        write(819+mu,*) gamma
        gamma = 3.1415926 * gamma / 2.d0
        
        write(6,*) 'gamma*pi/2=',gamma*pi/2
        !
        ! the factor 2 comes from the factor sqrt(hbar/2/M/omega) that appears
        ! in the definition of the electron-phonon matrix element g
        ! The sqrt(1/M) factor is actually hidden into the normal modes
        !
        ! gamma = \pi \sum_k\sum_{i,j} \delta(e_{k,i}-Ef) \delta(e_{k+q,j}-Ef)
        !         | \sum_mu z(mu,nu) <psi_{k+q,j}|dvscf_q(mu)*psi_{k,i}> |^2
        ! where z(mu,nu) is the mu component of normal mode nu (z = dyn)
        ! gamma(nu) is the phonon linewidth of mode nu
        !
        ! The factor N(Ef)^2 that appears in most formulations of el-ph interact
        ! is absent because we sum, not average, over the Fermi surface.
        ! The factor 2 is provided by the sum over spins
        !
        
        IF (SQRT (ABS (w2 (nu) ) ) .GT.eps) THEN
           ! lambda is the adimensional el-ph coupling for mode nu:
           ! lambda(nu)= gamma(nu)/(pi N(Ef) \omega_{q,nu}^2)
           lambda = gamma / 3.1415926 / w2 (nu) / dosef
        ELSE
           lambda = 0.0
        ENDIF
        ! 3.289828x10^6 is the conversion factor from Ry to GHz
        WRITE (6, 9010) nu, lambda, gamma * 3.289828d6
        IF (iuelph.NE.0) WRITE (iuelph, 9010) nu, lambda, gamma * &
             3.289828d6
     ENDDO
  ENDDO
  

9000 FORMAT(5x,'Gaussian Broadening: ',f7.3,' Ry, ngauss=',i4)
9005 FORMAT(5x,'DOS =',f10.6,' states/spin/Ry/Unit Cell at Ef=', &
          &       f10.6,' eV')
9006 FORMAT(5x,'double delta at Ef =',f10.6)
9010 FORMAT(5x,'lambda(',i2,')=',f8.4,'   gamma=',f8.2,' GHz')
  !
  !
  IF (iuelph.NE.0) CLOSE (unit = iuelph)
  RETURN
  

     
     !          call star_q(x_q(1,iq), at, bg, nsym , s , invs , nq, sxq, &
     !               isq, imq, .FALSE. )
     
     
  return
END SUBROUTINE ELPHSUM_WANNIER
   
