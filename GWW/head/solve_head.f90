! FOR GWW
!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Taken from Phonon code
! Modified by P. Umari and G. Stenuit
!
!-----------------------------------------------------------------------
subroutine solve_head
  !
!#ifdef __GWW
  !-----------------------------------------------------------------------
  !
  !calculates the head and wings of the dielectric matrix
  !
  USE ions_base,             ONLY : nat
  USE io_global,             ONLY : stdout, ionode,ionode_id
  USE io_files,              ONLY : prefix, iunigk, find_free_unit, diropn
  use pwcom
  USE check_stop,            ONLY : max_seconds
  USE wavefunctions_module,  ONLY : evc
  USE kinds,                 ONLY : DP
  USE becmod,                ONLY : bec_type, becp
  USE uspp_param,            ONLY : nhm
  USE uspp,                  ONLY : nkb, vkb, okvan
  use phcom,                 ONLY : dpsi, dvpsi, reduce_io, fildrho, iudrho, lrdrho, lgamma, &
                                    nksq, lrwfc, iuwfc, npwq, nbnd_occ, lrebar, iuebar
  USE wannier_gw,            ONLY : n_gauss, omega_gauss, grid_type
  USE control_ph,            ONLY : tr2_ph
  USE realus,                ONLY : adduspos_r
  USE gvect,                 ONLY : ig_l2g, mill_g
  USE mp_global,             ONLY : inter_pool_comm, intra_pool_comm
  USE mp,                    ONLY : mp_sum, mp_barrier, mp_bcast
  USE becmod,                ONLY : calbec
  USE symme,                 ONLY : symmatrix, crys_to_cart
  USE fft_base,             ONLY : dffts, dfftp
  USE fft_interfaces,       ONLY : fwfft, invfft

  implicit none

  real(DP) ::  thresh, anorm, averlt, dr2
  ! thresh: convergence threshold
  ! anorm : the norm of the error
  ! averlt: average number of iterations
  ! dr2   : self-consistency error
  real(DP), allocatable :: h_diag (:,:), eprec(:)
  ! h_diag: diagonal part of the Hamiltonian
  ! eprec : array fo preconditioning


  complex(DP) , allocatable :: auxg (:), aux1 (:),  ps (:,:)

  complex(DP), EXTERNAL :: ZDOTC      ! the scalar product function

  logical :: conv_root, exst
  ! conv_root: true if linear system is converged

  integer :: kter, iter0, ipol,jpol, ibnd, jbnd, iter, lter, &
       ik, ig, irr, ir, is, nrec, ios
  ! counters
  integer :: ltaver, lintercall

  real(DP) :: tcpu, get_clock
  ! timing variables

  character (len=256) :: flmixdpot
  ! the name of the file with the mixing potential

  external ch_psi_all, cg_psi

  REAL(kind=DP), ALLOCATABLE :: head(:),head_tmp(:)
  COMPLEX(kind=DP) :: sca, sca2
  REAL(kind=DP), ALLOCATABLE :: x(:),w(:), freqs(:)
  COMPLEX(kind=DP), ALLOCATABLE :: e_head(:,:)!wing of symmetric dielectric matrix (for G of local processor)
  COMPLEX(kind=DP), ALLOCATABLE :: e_head_g(:),e_head_g_tmp(:,:)
  COMPLEX(kind=DP), ALLOCATABLE :: e_head_pol(:,:,:)
  INTEGER :: i,j,k, iun
  REAL(kind=DP) :: ww, weight
  COMPLEX(kind=DP), ALLOCATABLE :: tmp_g(:)
  COMPLEX(kind=DP), ALLOCATABLE :: psi_v(:,:), prod(:), becpd(:,:)
  COMPLEX(kind=DP), ALLOCATABLE :: pola_charge(:,:,:)
  COMPLEX(kind=DP), ALLOCATABLE :: dpsi_ipol(:,:,:)
  REAL(kind=DP), ALLOCATABLE :: epsilon_g(:,:,:)
  INTEGER :: i_start,idumm,idumm1,idumm2,idumm3,ii
  REAL(kind=DP) :: rdumm

  write(stdout,*) 'Routine solve_head'
  call flush_unit(stdout)


  allocate(e_head(ngm,n_gauss+1))
  allocate(e_head_pol(ngm,n_gauss+1,3))
  e_head(:,:) =(0.d0,0.d0)
  allocate(x(2*n_gauss+1),w(2*n_gauss+1), freqs(n_gauss+1))
  allocate(head(n_gauss+1),head_tmp(n_gauss+1))
  head(:)=0.d0
  allocate(psi_v(dffts%nnr, nbnd), prod(dfftp%nnr))
  allocate (becpd (nkb, nbnd), tmp_g(ngm))
  allocate( pola_charge(dfftp%nnr,nspin,3))
  allocate(dpsi_ipol(npwx,nbnd,3),epsilon_g(3,3,n_gauss+1))
  epsilon_g(:,:,:)=0.d0
  e_head_pol(:,:,:)=0.d0
!setup Gauss Legendre frequency grid
!IT'S OF CAPITAL IMPORTANCE TO NULLIFY THE FOLLOWING ARRAYS
  x(:)=0.d0
  w(:)=0.d0
  if(grid_type==0) then
     call legzo(n_gauss*2+1,x,w)
     freqs(1:n_gauss+1)=-x(n_gauss+1:2*n_gauss+1)*omega_gauss
  else
     call legzo(n_gauss,x,w)
     freqs(1) = 0.d0
     freqs(2:n_gauss+1)=(1.d0-x(1:n_gauss))*omega_gauss/2.d0
  endif
  do i=1,n_gauss+1
     write(stdout,*) 'Freq',i,freqs(i)
  enddo
  CALL flush_unit( stdout )

  deallocate(x,w)
  head(:)=0.d0

  if (lsda) call errore ('solve_head', ' LSDA not implemented', 1)

  call start_clock ('solve_head')

  allocate (auxg(npwx))
  allocate (aux1(dffts%nnr))
  allocate (ps  (nbnd,nbnd))
  ps (:,:) = (0.d0, 0.d0)
  allocate (h_diag(npwx, nbnd))
  allocate (eprec(nbnd))

  IF (ionode .AND. fildrho /= ' ') THEN
     INQUIRE (UNIT = iudrho, OPENED = exst)
     IF (exst) CLOSE (UNIT = iudrho, STATUS='keep')
     CALL diropn (iudrho, TRIM(fildrho)//'.E', lrdrho, exst)
  end if
  !
  ! if q=0 for a metal: allocate and compute local DOS at Ef
  !
  if (degauss.ne.0.d0.or..not.lgamma) call errore ('solve_e', &
       'called in the wrong case', 1)
  !
  if (reduce_io) then
     flmixdpot = ' '
  else
     flmixdpot = 'mixd'
  endif
  !
  !   only one iteration is required
  !
  !   rebuild global miller index array
  !
  allocate(mill_g(3,ngm_g))
  mill_g(:,:) = 0
  do ig = 1, ngm
     i = nint( g(1,ig)*bg(1,1) + g(2,ig)*bg(2,1) + g(3,ig)*bg(3,1) )
     j = nint( g(1,ig)*bg(1,2) + g(2,ig)*bg(2,2) + g(3,ig)*bg(3,2) )
     k = nint( g(1,ig)*bg(1,3) + g(2,ig)*bg(2,3) + g(3,ig)*bg(3,3) )
     mill_g (1, ig_l2g(ig)) = i
     mill_g (2, ig_l2g(ig)) = j
     mill_g (3, ig_l2g(ig)) = k
  end do
  call mp_sum( mill_g )
  !
  if(ionode) then

     inquire(file=trim(prefix)//'.head_status', exist = exst)
     if(.not. exst) then
        i_start=1
     else
        iun =  find_free_unit()
        open( unit= iun, file=trim(prefix)//'.head_status', status='old')
        read(iun,*) i_start
        close(iun)
        if(i_start<1 .or. i_start>=(n_gauss+1)) then
           i_start=1
        else
           i_start=i_start+1
        endif
     endif
  endif
  call mp_bcast(i_start,ionode_id)

  do i=i_start,n_gauss+1
     write(stdout,*) 'Freq',i
     call flush_unit(stdout)
     pola_charge(:,:,:)=(0.d0,0.d0)
     if (nksq.gt.1) rewind (unit = iunigk)
     do ik = 1, nksq
        !
        write(stdout,*) 'ik:', ik
        call flush_unit(stdout)
        !
        weight = wk (ik)
        ww = fpi * weight / omega

        if (lsda) current_spin = isk (ik)
        if (nksq.gt.1) then
           read (iunigk, err = 100, iostat = ios) npw, igk
100        call errore ('solve_e', 'reading igk', abs (ios) )
        endif
        !
        ! reads unperturbed wavefuctions psi_k in G_space, for all bands
        !
        if (nksq.gt.1) call davcio (evc, lrwfc, iuwfc, ik, - 1)
        npwq = npw
        call init_us_2 (npw, igk, xk (1, ik), vkb)

        !trasnform valence wavefunctions to real space
        do ibnd=1,nbnd
           psi_v(:,ibnd) = ( 0.D0, 0.D0 )
           psi_v(nls(igk(1:npw)),ibnd) = evc(1:npw,ibnd)
           CALL invfft ('Wave', psi_v(:,ibnd), dffts)
        enddo !do ibnd=1,nbnd
        !
        ! compute the kinetic energy
        !
        do ig = 1, npwq
           g2kin (ig) = ( (xk (1,ik ) + g (1,igk (ig)) ) **2 + &
                (xk (2,ik ) + g (2,igk (ig)) ) **2 + &
                (xk (3,ik ) + g (3,igk (ig)) ) **2 ) * tpiba2
        enddo !do ig = 1, npwq
        !
        dpsi_ipol(:,:,:)=(0.d0,0.d0)
        do ipol = 1,3
           write(stdout,*) 'ipol:', ipol
           call flush_unit(stdout)
           !
           ! computes/reads P_c^+ x psi_kpoint into dvpsi array
           !
           call dvpsi_e (ik, ipol)
           !
           ! Orthogonalize dvpsi to valence states: ps = <evc|dvpsi>
           !
           CALL ZGEMM( 'C', 'N', nbnd_occ (ik), nbnd_occ (ik), npw, &
                (1.d0,0.d0), evc(1,1), npwx, dvpsi(1,1), npwx, (0.d0,0.d0), &
                ps(1,1), nbnd )
#ifdef __MPI
           call mp_sum(ps)
           !!!call reduce (2 * nbnd * nbnd_occ (ik), ps)
#endif
           ! dpsi is used as work space to store S|evc>
           !
           !!!CALL ccalbec (nkb, npwx, npw, nbnd_occ(ik), becp, vkb, evc)
           call calbec(npw, vkb, evc, becp, nbnd_occ(ik))
           CALL s_psi (npwx, npw, nbnd_occ(ik), evc, dpsi)
           !
           ! |dvpsi> = - (|dvpsi> - S|evc><evc|dvpsi>)
           ! note the change of sign!
           !
           CALL ZGEMM( 'N', 'N', npw, nbnd_occ(ik), nbnd_occ(ik), &
                (1.d0,0.d0), dpsi(1,1), npwx, ps(1,1), nbnd, (-1.d0,0.d0), &
                dvpsi(1,1), npwx )
           ! starting threshold for the iterative solution of the linear
           thresh = tr2_ph!ATTENZIONE
           !
           ! iterative solution of the linear system (H-e)*dpsi=dvpsi
           ! dvpsi=-P_c+ (dvbare)*psi
           !
           do ibnd = 1, nbnd_occ (ik)
              do ig = 1, npw
                 auxg (ig) = g2kin (ig) * evc (ig, ibnd)
              enddo
              eprec (ibnd) = 1.35d0*ZDOTC(npwq,evc(1,ibnd),1,auxg,1)
           enddo ! do ibnd = 1, nbnd_occ (ik)
#ifdef __MPI
           call mp_sum(eprec)
           !!!call reduce (nbnd_occ (ik), eprec)
#endif
           do ibnd = 1, nbnd_occ (ik)
              do ig = 1, npw
                 h_diag(ig,ibnd)=(1.d0/max(1.0d0,g2kin(ig)/eprec(ibnd)))**2.d0
              enddo
           enddo !do ibnd = 1, nbnd_occ (ik)
           !
           conv_root = .true.
           !TD put also imaginary frequency +iw
           call cgsolve_all_imfreq (ch_psi_all,cg_psi,et(1,ik),dvpsi,dpsi, &
                h_diag,npwx,npw,thresh,ik,lter,conv_root,anorm,nbnd_occ(ik),freqs(i))
           !dvpsi is NOT currupted on exit
           if (.not.conv_root) WRITE( stdout, "(5x,'kpoint',i4,' ibnd',i4, &
                &         ' solve_head: root not converged ',e10.3)") ik &
                &, ibnd, anorm
           !
           dpsi_ipol(:,:,ipol)=dpsi(:,:)
           !it calculates the wings \epsilon(G=0,G/=0)
           ! = <dpsi|exp(iGx)|psi_v>
           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           !if required calculates betae for dpsi
           if (okvan) call calbec(npw, vkb, dpsi, becpd, nbnd_occ(ik))
           !!!if(okvan) CALL ccalbec (nkb, npwx, npw, nbnd_occ(ik), becpd, vkb, dpsi)
           ! cycle on bands
           do ibnd=1,nbnd
              !      fft trasform dpsi to real space
              prod(:) = ( 0.D0, 0.D0 )
              prod(nls(igk(1:npw))) = dpsi(1:npw,ibnd)
              !
              CALL invfft ('Wave', prod, dffts)
              !      product dpsi * psi_v
              prod(1:dffts%nnr)=conjg(prod(1:dffts%nnr))*psi_v(1:dffts%nnr,ibnd)
              if(doublegrid) then
                 call cinterpolate(prod,prod,1)
              endif
              !      add us part if required
              if(okvan) call adduspos_r(prod,becpd(:,ibnd),becp%k(:,ibnd))
              !
              pola_charge(:,1,ipol)=pola_charge(:,1,ipol)+prod(:)*ww
              !
           enddo ! do ibnd=1,nbnd
           !
        enddo ! do ipol = 1,3 (on polarization) line 224
        !
        ! update epsilon tensor
        !
        do ipol = 1, 3
           nrec = (ipol - 1) * nksq + ik
           call davcio (dvpsi, lrebar, iuebar, nrec, - 1)
           do jpol = 1, 3
              !
              do ibnd = 1, nbnd_occ (ik)
                 !
                 !  this is the real part of <DeltaV*psi(E)|DeltaPsi(E)>
                 !
                 epsilon_g(ipol,jpol,i)=epsilon_g(ipol,jpol,i)-4.d0*ww* DBLE( &
                      ZDOTC (npw, dvpsi (1, ibnd), 1, dpsi_ipol (1, ibnd,ipol), 1) )
              enddo ! do ibnd = 1, nbnd_occ (ik)
           enddo ! do jpol = 1, 3
        enddo ! do ipol = 1, 3
        !
        !
     enddo ! do ik = 1, nksq (on k-points) line 188
     !
     !
     !      fft trasform to g space
     !      extract terms
#ifdef __MPI
     call mp_sum ( pola_charge, inter_pool_comm )
     !!!call poolreduce (2 * 3 * dfftp%nnr *nspin, pola_charge)
     call psyme (pola_charge)
#else
     call syme (pola_charge)
#endif
     !
     do ipol=1,3
        CALL fwfft ('Dense', pola_charge(:,1,ipol), dfftp)
        !
        tmp_g(:)=(0.d0,0.d0)
        tmp_g(gstart:ngm)=pola_charge(nl(gstart:ngm),1,ipol)
        !
        sca=(0.d0,0.d0)
        do ig=1,ngm
           sca=sca+conjg(tmp_g(ig))*tmp_g(ig)
        enddo
        call mp_sum(sca)
        write(stdout,*) 'POLA SCA', sca
        ! loop on frequency
        do ig=gstart,ngm
           e_head_pol(ig,i,ipol)=e_head_pol(ig,i,ipol)-4.d0*tmp_g(ig)
        enddo
     enddo ! do ipol=1,3
     !
#ifdef __MPI
     call mp_sum (epsilon_g(:,:,i), intra_pool_comm)
     !!!call reduce (9, epsilon_g(:,:,i))
     call mp_sum ( epsilon_g(:,:,i), inter_pool_comm )
     !!!call poolreduce (9, epsilon_g(:,:,i))
#endif
     !
     !   symmetrize
     !
     WRITE( stdout,'(/,10x,"Unsymmetrized in crystal axis ",/)')
     WRITE( stdout,'(10x,"(",3f15.5," )")') ((epsilon_g(ipol,jpol,i),&
     &                                ipol=1,3),jpol=1,3)
     !
     call crys_to_cart (epsilon_g(:,:,i))
     call symmatrix (epsilon_g(:,:,i) )
     !
     !    pass to cartesian axis
     !
     WRITE( stdout,'(/,10x,"Symmetrized in cartesian axis ",/)')
     WRITE( stdout,'(10x,"(",3f15.5," )")') ((epsilon_g(ipol,jpol,i),&
     &                                ipol=1,3),jpol=1,3)
     !
     ! add the diagonal part
     !  and print the result
     !
     WRITE( stdout, '(/,10x,"Dielectric constant in cartesian axis ",/)')

     WRITE( stdout, '(10x,"(",3f18.9," )")') ((epsilon_g(ipol,jpol,i), ipol=1,3), jpol=1,3)

     ! reminder i runs over freq (do i=i_start,n_gauss+1, line 183)
     head(i) = (epsilon_g(1,1,i)+epsilon_g(2,2,i)+epsilon_g(3,3,i))/3.d0
     e_head(:,i)=(e_head_pol(:,i,1)+e_head_pol(:,i,2)+e_head_pol(:,i,3))/3.0
     !
     ! TD writes on files
     write(stdout,*) 'ionode=', ionode
     call flush_unit(stdout)
     if(ionode) then
        !
        write(stdout,*) 'HEAD:',freqs(i),head(i)
        !
        write(stdout,*) 'E_HEAD :', i, e_head(2, i)
        write(stdout,*) i,e_head_pol(2,i,1)
        write(stdout,*) i,e_head_pol(2,i,2)
        write(stdout,*) i,e_head_pol(2,i,3)
        !
     endif
     !
     call flush_unit(stdout)
     !
     ! writes on file
     !
     if(ionode) then
        iun =  find_free_unit()
        if(i==1) then
           open( unit= iun, file=trim(prefix)//'.head', status='unknown',form='unformatted')
           write(iun) n_gauss
           write(iun) omega_gauss
           write(iun) freqs(1:n_gauss+1)
           write(iun) head(1:n_gauss+1)
        else
           open( unit= iun, file=trim(prefix)//'.head', status='old',position='rewind',form='unformatted')
           read(iun) idumm
           read(iun) rdumm
           read(iun) head_tmp(1:n_gauss+1)
           read(iun) head_tmp(1:n_gauss+1)
           head(1:i-1)=head_tmp(1:i-1)
           rewind(iun)
           write(iun) n_gauss
           write(iun) omega_gauss
           write(iun) freqs(1:n_gauss+1)
           write(iun) head(1:n_gauss+1)
        endif
        close(iun)
     endif
     !
     !write(stdout,*) 'file .head updated'
     !call flush_unit(stdout)
     !
     ! collect data
     allocate(e_head_g(ngm_g))
     e_head_g(:)=(0.d0,0.d0)
     !write(stdout,*) 'check1'
     !call flush_unit(stdout)
     !
     do ig = 1, ngm
        e_head_g( ig_l2g( ig ) ) = e_head_pol(ig,i,1)
     end do
     !write(stdout,*) 'check2'
     !call flush_unit(stdout)
     call mp_sum( e_head_g )
     !
     write(stdout,*) 'check3'
     call flush_unit(stdout)
     !
     write(stdout,*) 'check3'
     if(ionode) then
        iun =  find_free_unit()
        if(i==1) then
           open( unit= iun, file=trim(prefix)//'.e_head', status='unknown',form='unformatted')

           write(stdout,*) 'n_gauss=',n_gauss
           write(stdout,*) 'omega_gauss=',omega_gauss
           write(stdout,*) 'freqs(1:n_gauss+1)=', freqs(1:n_gauss+1)
           write(stdout,*) 'ngm_g=',ngm_g
           call flush_unit(stdout)
           write(iun) n_gauss
           write(iun) omega_gauss
           write(iun) freqs(1:n_gauss+1)
           write(iun) ngm_g
           write(stdout,*) 'size of mill_g=', size(mill_g)
           write(stdout,*) 'should be equal to 3*ngm_g=', 3*ngm_g
           write(stdout,*) 'mill_g(1,1)=', mill_g(1,1)
           write(stdout,*) 'mill_g(1,ngm_g)=', mill_g(1,ngm_g)
           write(stdout,*) 'mill_g(2,1)=', mill_g(2,1)
           write(stdout,*) 'mill_g(2,ngm_g)=', mill_g(2,ngm_g)
           write(stdout,*) 'mill_g(3,1)=', mill_g(3,1)
           write(stdout,*) 'mill_g(3,ngm_g)=', mill_g(3,ngm_g)
           call flush_unit(stdout)
           do ig=1,ngm_g
              write(iun) mill_g(1,ig), mill_g(2,ig), mill_g(3,ig)
           enddo
           write(stdout,*) 'coucou10'
           call flush_unit(stdout)
           do ig=1,ngm_g
              write(iun) e_head_g(ig)
           enddo
           write(stdout,*) 'check4 inside if(ionode) and if(i==1)'
           call flush_unit(stdout)
        else
           if(i/=(n_gauss+1)) then
              open( unit= iun, file=trim(prefix)//'.e_head', status='old',position='append',form='unformatted')
              do ig=1,ngm_g
                 write(iun) e_head_g(ig)
              enddo
           else ! si if i=n_gauss+1 (last freq)
              open( unit= iun, file=trim(prefix)//'.e_head', status='old',position='rewind',form='unformatted')
              read(iun) idumm
              read(iun) rdumm
              read(iun) head_tmp(1:n_gauss+1)
              read(iun) idumm
              do ig=1,ngm_g
                 read(iun) idumm1,idumm2,idumm3
              enddo
              allocate(e_head_g_tmp(n_gauss+1,ngm_g))
              do ii=1,n_gauss
                 do ig=1,ngm_g
                    read(iun) e_head_g_tmp(ii,ig)
                 enddo
              enddo
              e_head_g_tmp(n_gauss+1,:)=e_head_g(:)
              rewind(iun)
              write(iun) n_gauss
              write(iun) omega_gauss
              write(iun) freqs(1:n_gauss+1)
              write(iun) ngm_g
              do ig=1,ngm_g
                 write(iun) mill_g(1,ig), mill_g(2,ig), mill_g(3,ig)
              enddo
              do ig=1,ngm_g
                 write(iun) e_head_g_tmp(1:n_gauss+1,ig)
              enddo
              deallocate(e_head_g_tmp)
           endif ! if(i/=(n_gauss+1)) line 460
           !
        endif ! if(i==1) then line 447
        write(stdout,*) 'coucou'
        call flush_unit(stdout)
        close(iun)
        !
     endif ! if(ionode) then line 445
     !
     write(stdout,*) 'Before deallocate(e_head_g)'
     call flush_unit(stdout)
     !
     deallocate(e_head_g)
     !
     write(stdout,*) 'file .e_head updated'
     call flush_unit(stdout)
     !
     if(ionode) then
        !
        ! write which Freq has been done in .head_status
        ! so if it crashes, it will restart from this Freq.+1
        !
        iun =  find_free_unit()
        open( unit= iun, file=trim(prefix)//'.head_status', status='unknown')
        write(iun,*) i
        close(iun)
        !
     endif ! if(ionode) then
     !
     !write(stdout,*) 'file .head_status updated'
     !call flush_unit(stdout)
     !write(stdout,*) 'before mp_barrier'
     !call flush_unit(stdout)
     !call mp_barrier
     !write(stdout,*) 'mp_barrier done'
     !write(stdout,*) 'Freq=', i, ' DONE'
     !call flush_unit(stdout)
     !
  enddo ! do i=i_start,n_gauss+1 (on the Freq) line 183
  !
  deallocate (mill_g)
  deallocate (eprec)
  deallocate (h_diag)
  deallocate (ps)
  deallocate (aux1)
  deallocate (auxg)
  deallocate(psi_v, prod, becpd)
  deallocate(pola_charge)

  deallocate(head,freqs)
  deallocate(e_head, tmp_g)
  deallocate(dpsi_ipol,epsilon_g)
  deallocate(e_head_pol)


  call stop_clock ('solve_head')
  !
!#endif __GWW
  !
  return
end subroutine solve_head
