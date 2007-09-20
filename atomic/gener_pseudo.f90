!
! Copyright (C) 2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!--------------------------------------------------------------------------
subroutine gener_pseudo
  !--------------------------------------------------------------------------
  !
  !     This routine generate a pseudopotential in separable form
  !     It can be of NC type or of US type
  !     Multiple projections are allowed.
  !     Spin-orbit split pseudopotentials are also available.
  !     NB: bmat indices are as in the Vanderbilt paper PRB (1990)
  !
  !     The output of the routine are:
  !
  !     phis: the pseudo wavefunctions
  !     betas: the nonlocal projectors
  !     bmat:  the pseudopotential coefficients
  !     qq:    the integrals of the q functions
  !     qvan:  the augmentation functions
  !     vpsloc: the local pseudopotential
  !     chis:   auxiliary functions
  !
  !
  !     The construction of a PAW dataset can also be done (experimental)
  !      
  use kinds, only: dp
  use radial_grids, only: ndmx
  use ld1_parameters, only: nwfsx
  use io_global, only : stdout, ionode, ionode_id
  use mp,        only : mp_bcast
  use ld1inc, only: grid, lls, jjs, new, psipaw, psi, ecutwfc, ecutrho, &
                    psipsus, tm, ocs, phis, els, nwfs, nspin, rel, nlcc, &
                    file_chi, file_beta, file_wfcncgen, file_qvan, &
                    file_wfcusgen, file_wfcaegen, psccharge, aeccharge, &
                    qvan, qq, bmat, ddd, betas, nbeta, ikk, pseudotype, &
                    pawsetup, zval, vpsloc, vpot, vnl, lpaw, rcloc, rcutus, &
                    enl, enls, rcut, chis, nstoae, paw_rmatch_augfun,&
                    lnc2paw, rcutnc2paw, rhos, which_paw_augfun
  use atomic_paw, only : us2paw, paw2us
  implicit none

  integer ::   &
       ik,    &  ! the point corresponding to rc
       ikus,  &  ! the point corresponding to rc ultrasoft
       ikloc, &  ! the point corresponding to rc local
       ns,    &  ! counter on pseudo functions
       ns1,   &  ! counter on pseudo functions
       ib,jb, &  ! counter on beta functions
       nnode, &  ! the number of nodes of phi
       lam       ! the angular momentum

  real(DP) ::    &
       xc(8),        &  ! parameters of bessel functions
       psi_in(ndmx),  &  ! the all_electron wavefunction
       gi(ndmx,2),    &  ! auxiliary to compute the integrals
       occ,          &
       sum, db, work(nwfsx) ! work space

  real(DP), allocatable :: &
       b(:,:), binv(:,:) ! the B matrix and its inverse

  real(DP) ::    &
       aekin(nwfsx,nwfsx),  & ! AE kinetic energies
       pskin(nwfsx,nwfsx),  & ! PS kinetic energies
       kindiff(nwfsx,nwfsx)   ! AE-PS k.e.

  real(DP), external ::    &
       int_0_inf_dr    ! the function calculating the integral 

  integer :: &
       m, n, l, n1, n2, nwf0, nst, ikl, imax, iwork(nwfsx), &
       is, nbf, nc, ios, ind

  character(len=5) :: indqvan

  logical :: &
       lbes4     ! use 4 Bessel functions expansion

  ! additional vars for paw
  real(DP) :: vpotpaw (ndmx) ! total potential to be used for PAW 
                             ! generation (normally the AE potential)
  integer  :: iknc2paw       ! point in rgrid closer to rcutnc2paw

  if (lpaw) then
     write(stdout, &
          '(/,5x,21(''-''),'' Generating PAW atomic setup '',20(''-''),/)')
  elseif (pseudotype == 1.or.pseudotype == 2) then
     write(stdout, &
          '(/,5x,21(''-''),'' Generating NC pseudopotential '',21(''-''),/)')
  elseif (pseudotype == 3) then
     write(stdout, &
          '(/,5x,21(''-''),'' Generating US pseudopotential '',21(''-''),/)')
  else
     call errore('gener_pseudo','pseudotype not programmed',1)
  endif
  if (pseudotype == 1.and.rel == 2) call errore('gener_pseudo', &
       'not programmed' ,2)
  if (pseudotype /= 3.and. lpaw) call errore('gener_pseudo', &
       'please start from a US for generating a PAW dataset' ,pseudotype)
  psipaw=0.0_dp
  phis=0.0_dp
  !
  !   compute the local potential from the all-electron potential
  !
  call pseudovloc ( )
  !
  !   initialize total potential for PAW generation
  if (lpaw) then
     if (.not.lnc2paw) then
        vpotpaw(1:grid%mesh) = vpot(1:grid%mesh,1)
     else
        vpotpaw(1:grid%mesh) = vpsloc(1:grid%mesh)
     end if
  endif
  !
  !   if nlcc is true compute here the core charge
  !   the core charge is needed also for the PAW dataset
  !
  if (nlcc .or. lpaw) call set_rho_core
  !
  !   set the appropriate energies and the correspondence all-electron
  !   pseudo
  !
  do n=1,nwfs
     if (enls(n) == 0.0_dp) enls(n)=enl(nstoae(n))
  enddo
  !
  ! Set the all-electron wavefunctions, calculating those at user supplied
  ! energies. The wavefunctions are written on file at this point, so
  ! the user can check them also when the pseudopotential generation is  
  ! unsuccessful
  ! 
  do ns=1,nbeta
     ik=0
     nwf0=nstoae(ns)
     do n=1,grid%mesh
        if (grid%r(n).lt.rcut(ns)) ik=n
     enddo
     if (mod(ik,2) == 0) ik=ik+1
     if (new(ns)) then
        call set_psi_in(ik,lls(ns),jjs(ns),enls(ns),psipaw(1,ns))
     else
        psipaw(:,ns)=psi(:,1,nwf0)
     endif
  enddo

  if (file_wfcaegen .ne. ' ') then
     if (ionode) &
        open(unit=19,file=file_wfcaegen, status='unknown', iostat=ios, err=800)
800  call mp_bcast(ios, ionode_id)
     call errore('gener_pseudo','opening file '//file_wfcaegen,abs(ios))
     if (ionode) then
        do n=1,grid%mesh
           write(19,'(15f12.6)') grid%r(n), (psipaw(n,ns), ns=1,nwfs)
        enddo
        close(19)
     endif
  endif
  !
  !   compute the pseudowavefunctions by expansion in spherical
  !   bessel function before r_c
  !
  ecutrho=0.0_dp
  ecutwfc=0.0_dp
  do ns=1,nbeta
     lam=lls(ns)
     nst=(lam+1)*2
     nwf0=nstoae(ns)
     !    
     !  compute the ik closer to r_cut, r_cutus, rcloc
     !
     ik=0
     ikus=0
     ikloc=0
     do n=1,grid%mesh
        if (grid%r(n).lt.rcut(ns)) ik=n
        if (grid%r(n).lt.rcutus(ns)) ikus=n
        if (grid%r(n).lt.rcloc) ikloc=n
     enddo
     if (mod(ik,2) == 0) ik=ik+1
     if (mod(ikus,2) == 0) ikus=ikus+1
     if (mod(ikloc,2) == 0) ikloc=ikloc+1
     if (lnc2paw) then
        do n=1,grid%mesh
           if (grid%r(n).lt.rcutnc2paw(ns)) iknc2paw=n
        end do
        if (mod(iknc2paw,2) == 0) iknc2paw=iknc2paw+1
     end if
     if (ikus.gt.grid%mesh) call errore('gener_pseudo','ik is wrong ',1)
     if (pseudotype == 3) then
        ikk(ns)=max(ikus+10,ikloc+5)
     else
        ikk(ns)=max(ik+10,ikloc+5)
     endif

     if (new(ns)) then
        occ=1.0_DP
     else
        occ=ocs(ns)
     endif
     !
     !   save the all-electron function for the PAW setup
     !
     psi_in(1:grid%mesh) = psipaw(1:grid%mesh,ns) 
     !
     !  compute the phi functions
     !
     if (lpaw.and.lnc2paw) then
        ! first compute possibly harder NC pseudowfcts to be
        ! used as AE reference for PAW generation
        nnode=0
       !call compute_phi(lam,ik,      psi_in,phis(1,ns),xc,1,occ,enls(ns),els(ns)) CURRENT
        call compute_phi(lam,iknc2paw,psi_in,phis(1,ns),xc,1,occ,enls(ns),els(ns))
       !call compute_phi(lam,iknc2paw,ik,nwf0,ns,xc,1,nnode,ocs(ns)) PAW
       !call compute_phi(lam,ik,      ik,nwf0,ns,xc,1,nnode,ocs(ns)) OLD
        psipaw(1:grid%mesh,ns)=phis(1:grid%mesh,ns)
     endif
     !
     if (tm) then
        call compute_phi_tm(lam,ik,psi_in,phis(1,ns),1,xc,enls(ns),els(ns))
     else
        call compute_phi(lam,ik,psi_in,phis(1,ns),xc,1,occ,enls(ns),els(ns))
        ecutrho=max(ecutrho,8.0_dp*xc(6)**2)
     endif
     !
     !   US only on the components where ikus <> ik
     ! 
     psipsus(:,ns)=phis(:,ns) 
     if (ikus.ne.ik) then
        call compute_phius(lam,ikus,psipsus(1,ns),phis(1,ns),xc,1,els(ns))
        ecutwfc=max(ecutwfc,2.0_dp*xc(5)**2)
        lbes4=.true.
     else
        lbes4=.false.
        if (.not.tm) ecutwfc=max(ecutwfc,2.0_dp*xc(6)**2)
     endif
     if (tm.and.ik==ikus) then
        call compute_chi_tm(lam,ik,ikk(ns),phis(1,ns),chis(1,ns),xc,enls(ns))
     else
        call compute_chi(lam,ikk(ns),phis(1,ns),chis(1,ns),xc,enls(ns),lbes4)
     endif
  enddo


  !      do n=1,mesh
  !         write(stdout,'(5e15.7)') r(n),psipsus(n,1),chis(n,1),
  !     +                            psipsus(n,2),chis(n,2)
  !      enddo
  !      stop

  !
  !    for each angular momentum take the same integration point
  !
  do ns=1,nbeta
     do ns1=1,nbeta
        if (lls(ns) == lls(ns1).and.ikk(ns1).gt.ikk(ns)) &
             ikk(ns)=ikk(ns1)
     enddo
  enddo
  !
  !     construct B_{ij}
  !
  bmat=0.0_dp
  do ns=1,nbeta
     do ns1=1,nbeta
        if (lls(ns) == lls(ns1).and.abs(jjs(ns)-jjs(ns1)).lt.1.e-7_dp) then
           nst=(lls(ns)+1)*2
           ikl=ikk(ns1)
           do n=1,grid%mesh
              gi(n,1)=phis(n,ns)*chis(n,ns1)
           enddo
           bmat(ns,ns1)=int_0_inf_dr(gi,grid,ikl,nst)
        endif
     enddo
  enddo

  allocate ( b(nbeta, nbeta), binv(nbeta, nbeta) )

  if (pseudotype == 1) then
     !
     !     NC single-projector PP: construct the semilocal potential 
     !
     vnl=0.0_dp
     do ns=1,nbeta
        lam=lls(ns)
        if ( rel < 2 .or. lls(ns) == 0 .or. &
             abs(jjs(ns)-lls(ns)+0.5_dp) < 0.001_dp) then
           ind=1
        else if ( rel == 2 .and. lls(ns) > 0 .and. &
             abs(jjs(ns)-lls(ns)-0.5_dp) < 0.001_dp) then
           ind=2
        endif
        do n=1,ikk(ns)
           vnl(n,lam,ind) = chis(n,ns)/phis(n,ns)
        enddo
     enddo
     !
  else if (pseudotype == 2) then
     !
     !     symmetrize the B matrix
     !
     do ns=1,nbeta
        do ns1=1,ns-1
           bmat(ns,ns1)=0.5_dp*(bmat(ns,ns1)+bmat(ns1,ns))
           bmat(ns1,ns)=bmat(ns,ns1)
        enddo
     enddo
  end if
  !
  do ns=1,nbeta
     do ns1=1,nbeta
        b(ns,ns1)=bmat(ns,ns1)
     enddo
  enddo
  !
  !   compute the inverse of the matrix B_{ij}^-1
  !
  if (nbeta > 0) call invmat(nbeta, b, binv, db)
  !
  !   compute the beta functions
  !
  betas=0.0_dp
  do ns=1,nbeta
     do ns1=1,nbeta
        do n=1,grid%mesh
           betas(n,ns)=betas(n,ns)+ binv(ns1,ns)*chis(n,ns1)
        enddo
     enddo
  enddo
  deallocate (b, binv)
  !
  qq=0.0_dp
  if (pseudotype == 3) then
     !
     !    compute the Q functions
     !
     do ns=1,nbeta
        do ns1=1,ns
           ikl=max(ikk(ns),ikk(ns1))
           do n=1, ikl
              qvan(n,ns,ns1) = psipsus(n,ns) * psipsus(n,ns1) &
                   - phis(n,ns) * phis(n,ns1)
              gi(n,1)=qvan(n,ns,ns1)
           enddo
           do n=ikl+1,grid%mesh
              qvan(n,ns,ns1)=0.0_dp
           enddo
           !
           !     and puts its integral in qq
           !
           if (lls(ns) == lls(ns1).and.abs(jjs(ns)-jjs(ns1)).lt.1.e-8_dp) then
              nst=(lls(ns)+1)*2
              qq(ns,ns1)=int_0_inf_dr(gi,grid,ikk(ns),nst)
           endif
           !
           !     set the bmat with the eigenvalue part
           !
           bmat(ns,ns1)=bmat(ns,ns1)+enls(ns1)*qq(ns,ns1)
           !
           !    Use symmetry of the n,ns1 indeces to set qvan and qq and bmat
           !
           if (ns.ne.ns1) then
              do n=1,grid%mesh
                 qvan(n,ns1,ns)=qvan(n,ns,ns1)
              enddo
              qq(ns1,ns)=qq(ns,ns1)
              bmat(ns1,ns)=bmat(ns1,ns)+enls(ns)*qq(ns1,ns)
           endif
        enddo
     enddo
     write(stdout,'(/5x,'' The bmat matrix'')')
     do ns1=1,nbeta
        write(stdout,'(6f12.5)') (bmat(ns1,ns),ns=1,nbeta)
     enddo
     write(stdout,'(/5x,'' The qq matrix'')')
     do ns1=1,nbeta
        write(stdout,'(6f12.5)') (qq(ns1,ns),ns=1,nbeta)
     enddo
  endif

  do is=1,nspin
     ddd(:,:,is)=bmat(:,:)
  enddo
  !
  !    generate a PAW dataset if required
  !
  if (lpaw) then
     if (lnc2paw) write (stdout,'(/5x,''WARNING: __PAW_FROM_NC__'')')
     !
     !symbol=atom_name(nint(zed))
     !
     ! compute kinetic energy differences, using:
     ! AE:   T |psi> = (e - Vae) |psi>
     ! PS:   T |phi> = (e - Vps) |phi> - |chi>
     do ns=1,nbeta
        do ns1=1,ns
           if (lls(ns)==lls(ns1).and.jjs(ns)==jjs(ns1)) then
              ikl=max(ikk(ns),ikk(ns1))
              nst=2*(lls(ns)+1)
              do n=1,ikl
                 gi(n,1)=psipaw(n,ns)*(enls(ns1)-vpotpaw(n))*psipaw(n,ns1)
              end do
              aekin(ns,ns1)=int_0_inf_dr(gi(1:grid%mesh,1),grid,ikl,nst)
              do n=1,ikl
                 gi(n,1)=phis(n,ns)*( (enls(ns1)-vpsloc(n))*phis(n,ns1) - chis(n,ns1) )
              end do
              pskin(ns,ns1)=int_0_inf_dr(gi(1:grid%mesh,1),grid,ikl,nst)
           else
              aekin(ns,ns1)=0._dp
              pskin(ns,ns1)=0._dp
           end if
           kindiff(ns,ns1)=aekin(ns,ns1)-pskin(ns,ns1)
           kindiff(ns1,ns)=aekin(ns,ns1)-pskin(ns,ns1)
        end do
     end do
     !
     ! create the 'pawsetup' object containing the atomic setup for PAW
     call us2paw ( pawsetup,                                         &
          zval, grid, paw_rmatch_augfun, ikk,  &
          nbeta, lls, jjs, ocs, enls, els, rcutus, psipaw, phis, betas, &
          qvan, kindiff, nlcc, aeccharge, psccharge, vpotpaw, vpsloc, &
          which_paw_augfun)
     !
  endif
  !
  !    unscreen the local potential and the D coefficients
  !
  if (lpaw) then
     ! reread augmentation functions and descreened potentials from PAW
     call paw2us ( pawsetup, zval, grid, nbeta, lls, jjs, ikk, betas, &
                   qq, qvan, vpsloc, bmat, rhos, els, rcutus, pseudotype )
  else
     call descreening
  end if
  !
  !     print the main functions on files
  !
  if (file_beta .ne. ' ') then
     if (ionode) &
        open(unit=19,file=file_beta, status='unknown', iostat=ios, err=400)
400  call mp_bcast(ios, ionode_id)
     call errore('gener_pseudo','opening file '//file_beta,abs(ios))
     if (ionode) then
        do n=1,grid%mesh
           write(19,'(15f12.6)') grid%r(n), (betas(n,ns), ns=1,nbeta)
        enddo
        close(19)
     endif
  endif
  if (file_chi .ne. ' ') then
     if (ionode) &
        open(unit=19,file=file_chi, status='unknown', iostat=ios, err=600)
600  call mp_bcast(ios, ionode_id)
     call errore('gener_pseudo','opening file '//file_chi,abs(ios))
     if (ionode) then
        do n=1,grid%mesh
           write(19,'(15f12.6)') grid%r(n), (chis(n,ns), ns=1,nbeta)
        enddo
        close(19)
     endif
  endif
  if (file_qvan .ne. ' ') then
     do ns1=1,nbeta
        indqvan=' '
        if (ns1 < 10) then
           write(indqvan,'(".",i1)') ns1
        elseif (ns1 < 100) then
           write(indqvan,'(".",i2)') ns1
        else
           write(indqvan,'(".",i3)') ns1
        endif
        if (ionode) &
           open(unit=19,file=TRIM(file_qvan)//TRIM(indqvan), status='unknown', &
             iostat=ios, err=700)
700     call mp_bcast(ios, ionode_id)
        call errore('gener_pseudo','opening file '//file_qvan,abs(ios))
        if (ionode) then
           do n=1,grid%mesh
              write(19,'(15f12.6)') grid%r(n), (qvan(n,ns,ns1), ns=1,ns1)
           enddo
           close(19)
        endif
     enddo
  endif


  if (file_wfcncgen .ne. ' ') then
     if (ionode) & 
        open(unit=19,file=file_wfcncgen, status='unknown', iostat=ios,err=900)
900  call mp_bcast(ios, ionode_id)
     call errore('gener_pseudo','opening file '//file_wfcncgen,abs(ios))
     if (ionode) then
        do n=1,grid%mesh
           write(19,'(15f12.6)') grid%r(n), (psipsus(n,ns), ns=1,nwfs)
        enddo
        close(19)
     endif
  endif

  if (file_wfcusgen .ne. ' ') then
     if (ionode) & 
        open(unit=19,file=file_wfcusgen, status='unknown', iostat=ios,err=1000)
1000 call mp_bcast(ios, ionode_id)
     call errore('gener_pseudo','opening file '//file_wfcusgen,abs(ios))
     if (ionode) then
        do n=1,grid%mesh
           write(19,'(15f12.6)') grid%r(n), (phis(n,ns), ns=1,nwfs)
        enddo
        close(19)
     endif
  endif

  write(stdout,"(/,5x,19('-'),' End of pseudopotential generation ',19('-'),/)")

  return
end subroutine gener_pseudo
