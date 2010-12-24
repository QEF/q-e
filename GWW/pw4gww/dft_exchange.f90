! FOR GWW
! Author: P. Umari
! Modified by G. Stenuit
!
subroutine dft_exchange_k(nbnd_v,nbnd_s, ecutoff)
!this subroutine calculates the exchange
!energi term for every state and writes on disk
!if l_truncated_coulomb is false assumes that the
!first k-point is gamma and this will not  be used for the
!calculation
!IT WILL REQUIRE SOME WORK FOR IMPEMENTING SPIN

! #ifdef __GWW

  USE io_global,            ONLY : stdout, ionode, ionode_id
  USE io_files,             ONLY : find_free_unit, prefix, iunwfc, nwordwfc, iunigk
  USE kinds,                ONLY : DP
  USE cell_base,            ONLY: at, alat, tpiba, omega, tpiba2
  USE wannier_gw
  USE fft_base,             ONLY : dffts, dfftp
  USE fft_interfaces,       ONLY : fwfft, invfft
  USE gvect
  USE gvecs,              ONLY : nls, nlsm, doublegrid
  USE uspp
  USE wavefunctions_module, ONLY : psic, evc
  USE cell_base,            ONLY : at, bg, omega
  USE mp,                   ONLY : mp_sum, mp_bcast
  USE klist,                ONLY : nks, nkstot, wk, xk, nelec, nelup, neldw, ngk
  USE wvfct,                ONLY : nbnd, npwx, npw, igk, wg, et, g2kin, ecutwfc
  USE exx,                  ONLY : yukawa
  USE constants,            ONLY : e2, pi, tpi, fpi, rytoev
  USE realus,               ONLY : adduspos_r
  USE becmod,               ONLY : calbec
  USE buffers,              ONLY : save_buffer, get_buffer

  implicit none

   INTEGER, INTENT(in)  :: nbnd_v !number of valence states
   INTEGER, INTENT(in)  :: nbnd_s !number of states considered
   REAL(kind=DP), INTENT(in) :: ecutoff!cutoff in Rydberg for g sum


   COMPLEX (kind=DP), ALLOCATABLE :: wfcs_s(:,:), wfcs_vc(:,:)
   COMPLEX(DP), ALLOCATABLE :: bec_vc(:,:), bec_s(:,:)
   REAL(kind=DP), ALLOCATABLE ::fac(:)
   INTEGER :: ii, ik, ig, jj, i, it, iun, ir
   REAL(kind=DP) :: qsca,sca
   COMPLEX(kind=DP) :: csca
   REAL(kind=DP), ALLOCATABLE :: ex(:,:)
   COMPLEX(kind=DP), ALLOCATABLE :: prods(:),prod(:)
   REAL(kind=DP), ALLOCATABLE :: x(:),w(:), times(:)
   REAL(kind=DP) :: e_fermi
   LOGICAL, PARAMETER :: l_gauss = .true.!if true uses auxiliary function
   LOGICAL, PARAMETER :: l_cubic=.true.!uses formula for cubic cell
   REAL(kind=DP) ::alpha!decay of auxiliary function
   REAL(kind=DP), ALLOCATABLE :: f_gauss(:)
   REAL(kind=DP) :: f_gauss_int
   REAL(kind=DP) :: qx,qy,qz
   INTEGER :: ngm_max
   REAL(kind=DP) :: q(3),q1(3)


   allocate(wfcs_s(dffts%nnr,nbnd_s), wfcs_vc(dffts%nnr,nbnd))
   ALLOCATE( bec_vc( nkb, nbnd ) )
   ALLOCATE( bec_s( nkb, nbnd ) )
   allocate(fac(ngm))
   allocate(prods(dffts%nnr),prod(dfftp%nnr))
   allocate(ex(nbnd_s,2*n_gauss+2))!the last elements contains t=0 and conduction states
   allocate(f_gauss(nks))

!determine ngm_max
   ngm_max=0
   do ig=1,ngm
      if(gg(ig)*tpiba2 >= ecutoff) exit
      ngm_max=ngm_max+1
   enddo

   write(stdout,*) 'NGM MAX:', ngm_max, ngm



 !set up alpha
   alpha=5.d0

   ex(:,:)=0.d0

!calculate fermi level
   e_fermi=(minval(et(nbnd_v+1,1:nks))-maxval(et(nbnd_v,1:nks)))/2.d0 &
                    &  +maxval(et(nbnd_v,1:nks))

 write(stdout,*) 'E FERMI :', e_fermi
 call flush_unit(stdout)

!!!!
!!! TO BE REMOVED :
 !e_fermi=0.5
 !e_fermi=0.464403884076
 !!!!!
 !write(stdout,*) 'E FERMI from OLD VERSION:', e_fermi
 !call flush_unit(stdout)

!set up Gauss Legendre time grid
   allocate(x(2*n_gauss+1),w(2*n_gauss+1),times(2*n_gauss+1))
   x(:)=0.d0
   w(:)=0.d0
   if(grid_type==0) then
      call legzo(n_gauss*2+1,x,w)
      times(:)=-x(:)*tau_gauss
   else
      call legzo(n_gauss,x,w)
      times(n_gauss+1)=0.d0
      times(n_gauss+2:2*n_gauss+1)=(1.d0-x(1:n_gauss))*tau_gauss/2.d0
      times(1:n_gauss)=(-1.d0-x(1:n_gauss))*tau_gauss/2.d0
   endif
   do i=1,2*n_gauss+1
      write(stdout,*) 'TIME:', i,times(i)
   enddo

!sets factors terms
!sets factors terms
!this has already  been called   call exx_grid_init()


!k=Gamma which correspons to kpoint1
!set up
   !write(stdout,*) 'nks=', nks, ' and nkb=', nkb
   !call flush_unit(stdout)

   !!!! just for testing purpose, I set all igk to 0
   !!!!! igk(1:npwx)=0.0
   !!!!!!!!!!!!!!!!
   !!!write(stdout,*) 'ubound(igk(:))=', ubound(igk(:))
   !!!write(stdout,*) 'igk(:)=', igk(:)
   !!!call flush_unit(stdout)

   !!!IF ( nks > 1 ) REWIND( iunigk ) ! en theorie il devrait deja les connaitre
   !!npw = ngk(ik)

   !write(stdout,*) 'COUCOU0 and npw=', npw
   !call flush_unit(stdout)

   !!! the idea is to replace all the READ( iunigk ) npw, igk
   !!! by only READ( iunigk ) igk
   !!!IF ( nks > 1 ) READ( iunigk ) npw, igk
!!!! from PH, it seems we have to read in the old way !!!!
!!!   IF ( nks > 1 ) READ( iunigk ) npw, igk
   !!!IF ( nks > 1 ) READ( iunigk ) igk ! en theorie connu
   !
   CALL gk_sort(xk(1,1), ngm, g, ecutwfc/tpiba2, npw, igk, g2kin)
   !
   !write(stdout,*) 'igk(:)', igk(:)
   !call flush_unit(stdout)

   npw = ngk(1)

   !write(stdout,*) 'COUCOU1 npw=', npw
   !call flush_unit(stdout)

   IF ( nkb > 0 ) CALL init_us_2( npw, igk, xk(1,1), vkb )

   !write(stdout,*) 'COUCOU2'
   !call flush_unit(stdout)

!read in wavefunctions
   !!!! I have removed the comment, since I am not sure
   !!! the evc are defined already for each ik =>
   IF( nks > 1 ) CALL davcio( evc, nwordwfc, iunwfc, 1, -1 )
!   IF ( nks > 1 ) CALL get_buffer( evc, nwordwfc, iunwfc, 1)

!!!!
   !write(stdout,*) 'npw=', npw, ' nbnd_s=', nbnd_s
   !write(stdout,*)"ubound(1), ubound(2) of evc : ", ubound(evc(:,:),1), ubound(evc(:,:),2)
   !write(stdout,*) 'evc(1,1)=', evc(1,1) , ' and evc(npw,1)=', evc(npw,1)
   !write(stdout,*) 'evc(1,2)=', evc(1,2) , ' and evc(npw,2)=', evc(npw,2)
   !write(stdout,*) 'ubound(igk) : ', ubound(igk(:))
   !write(stdout,*) 'igk(1) = ', igk(1)
   !call flush_unit(stdout)
!writes the wavefunctions in real space coarse grid

   do ii=1,nbnd_s
      wfcs_s(:,ii) = ( 0.D0, 0.D0 )
      wfcs_s(nls(igk(1:npw)),ii) = evc(1:npw,ii)
      CALL invfft ('Wave', wfcs_s(:,ii), dffts)
   enddo
   IF ( nkb > 0 ) &
        CALL calbec(npw, vkb, evc, bec_s, nbnd)


!!!! loop on k-points
!!!   IF ( nks > 1 ) REWIND( iunigk ) ! probablement inutile
!loop on k-points valence

   !write(stdout,*) 'coucou'
   !call flush_unit(stdout)

!if required calculates the auxiliary functions
   if(l_gauss .and. .not.l_truncated_coulomb) then
      f_gauss(1)=0.d0
      if(.not.l_cubic) then
         do ik=2,nks
            f_gauss(ik)=0.d0
            do ig=1,ngm_max
               qsca = (-xk(1,1)+xk(1,ik)+g(1,ig))**2.d0 + &
                    (-xk(2,1)+xk(2,ik)+g(2,ig))**2.d0 + &
                    (-xk(3,1)+xk(3,ik)+g(3,ig))**2.d0
               qsca=qsca*tpiba2
               f_gauss(ik)=f_gauss(ik)+exp(-alpha*qsca)/qsca
            enddo
            call mp_sum(f_gauss(ik))
         enddo
         f_gauss(:)=f_gauss(:)*e2*fpi/omega
         f_gauss_int=2.d0*alpha*(pi/alpha)**(3.d0/2.d0)*e2*fpi/((2.d0*pi)**3.d0)
      else
         do ik=2,nks
            qx=-xk(1,1)+xk(1,ik)
            qy=-xk(2,1)+xk(2,ik)
            qz=-xk(3,1)+xk(3,ik)
            f_gauss(ik)=alat**2.d0/(3.d0-cos(tpi*qx)-cos(tpi*qy)-cos(tpi*qz))/2.d0
         enddo
         f_gauss_int=9.9774204d0/(tpiba**2.d0)
         f_gauss(:)=(fpi*e2/omega)*f_gauss(:)
         f_gauss_int=(fpi*e2/omega)*f_gauss_int
      endif
   endif

   !write(stdout,*) 'coucou2'
   !call flush_unit(stdout)

!---------------------------------------------------------
   do ik=1,nks
     write(stdout,*) 'IK :',ik
     call flush_unit(stdout)

      if(l_gauss .and. .not.l_truncated_coulomb) write(stdout,*) 'FGAUSS', ik, f_gauss(ik),f_gauss_int

      !!!! it seems we have to read igk in the old way !!!
      !!!!IF(nks > 1) READ( iunigk ) npw, igk
      !!! IF(nks > 1) READ( iunigk ) igk  ! en theorie connu
      !
      npw = ngk(ik)
      !
      CALL gk_sort(xk(1,ik), ngm, g, ecutwfc/tpiba2, npw, igk, g2kin)
      !
      !write(stdout,*) 'ik=', ik, ' and npw=', npw
      !write(stdout,*) 'igk= ', igk(:)
      !call flush_unit(stdout)

      IF ( nkb > 0 ) CALL init_us_2( npw, igk, xk(1,ik), vkb ) !ATTENZIONE
      IF ( nks > 1 ) CALL davcio( evc, nwordwfc, iunwfc, ik, -1 )
!      IF ( nks > 1 ) CALL get_buffer( evc, nwordwfc, iunwfc, ik)

      !write(stdout,*) 'ubound(1), ubound(2) of evc : ', ubound(evc(:,:),1), ubound(evc(:,:),2)
      !write(stdout,*) 'evc(1,1)=', evc(1,1) , ' and evc(npw,1)=', evc(npw,1)
      !write(stdout,*) 'evc(1,2)=', evc(1,2) , ' and evc(npw,2)=', evc(npw,2)
      !write(stdout,*) 'ubound(igk) : ', ubound(igk(:))
      !write(stdout,*) 'igk(1) = ', igk(1)
      !write(stdout,*) 'xk(1,ik) = ', xk(1,ik)
      !call flush_unit(stdout)

      do ii=1,nbnd
         wfcs_vc(:,ii) = ( 0.D0, 0.D0 )
         wfcs_vc(nls(igk(1:npw)),ii) = evc(1:npw,ii)
         CALL invfft ('Wave', wfcs_vc(:,ii), dffts)
      enddo
      IF ( nkb > 0 ) &
           CALL calbec(npw, vkb, evc, bec_vc, nbnd)
!calcolate q vector
      q(:)=xk(:,ik)-xk(:,1)
      q(:)=-q(:)*tpiba

      !write(stdout,*) 'q(:) = ', q(:)
      !call flush_unit(stdout)

!calcolate qsca
!calcolate fac(q)


      do ig=1,ngm_max
         qsca = (-xk(1,1)+xk(1,ik)+g(1,ig))**2.d0 + &
              (-xk(2,1)+xk(2,ik)+g(2,ig))**2.d0 + &
              (-xk(3,1)+xk(3,ik)+g(3,ig))**2.d0
         if(.not.l_truncated_coulomb) then
            if (qsca > 1.d-8) then
               fac(ig)=e2*fpi/(tpiba2*qsca + yukawa )
            else
               fac(ig)= 0.d0!- exxdiv ! & ! or rather something else (see F.Gygi)
 !                         - e2*fpi   ! THIS ONLY APPLYS TO HYDROGEN
               if (yukawa .gt. 1.d-8 ) then
                  fac(ig) = fac(ig) + e2*fpi/(tpiba2*qsca + yukawa )
               end if
            end if
         else
            if (qsca > 1.d-8) then
               fac(ig)=(e2*fpi/(tpiba2*qsca))*(1.d0-dcos(dsqrt(qsca)*truncation_radius*tpiba))
            else
               fac(ig)=e2*fpi*(truncation_radius**2.d0/2.d0)
            endif
         endif
      enddo
      fac(:)=fac(:)/omega

      !write(stdout,*) 'fac(1) = ', fac(1)
      !call flush_unit(stdout)

!loop on s states

      do ii=1,nbnd_s
         write(stdout,*) 'II:', ii
         call flush_unit(stdout)

!loop on negative times

!loop on c states
         do jj=nbnd_v+1,nbnd



!calcolate charge
            call start_clock('k prod')

            prods(:)=conjg(wfcs_s(:,ii))*wfcs_vc(:,jj)
            !write(stdout,*) 'nrxxs=', dffts%nnr
            !write(stdout,*) 'prods(1)=', prods(1), 'prodsnrxxs)=', prods(dffts%nnr)
            !call flush_unit(stdout)

            call stop_clock('k prod')

            call start_clock('k inte')
            IF ( doublegrid ) THEN
               CALL cinterpolate( prod, prods, 1 )
            ELSE
               prod(:)=prods(:)
            ENDIF
            call stop_clock('k inte')
!ULTRASOFT STUFF TO BE ADDED ATTENZIONE

            call start_clock('k add')

            if(okvan) call adduspos_r(prod,bec_s(:,ii),bec_vc(:,jj))
            call stop_clock('k add')

            call start_clock('k fft')
            CALL fwfft ('Dense', prod, dfftp)
            call stop_clock('k fft')


!calculate energy
            call start_clock('k sum')
            sca=0.d0
            do ig = 1, ngm_max
               sca=sca+conjg(prod(nl(ig)))*prod(nl(ig))*fac(ig)
            enddo
            call mp_sum(sca)

            !write(stdout,*) 'sca = ', sca
            !call flush_unit(stdout)

            call stop_clock('k sum')
            call start_clock('k other')
            do it=1,n_gauss+1
               if(it/=n_gauss+1) then
                  ex(ii,it)=ex(ii,it)+sca*wk(ik)*exp((et(jj,ik)-e_fermi)*times(it))
               else
                  ex(ii,2*n_gauss+2)=ex(ii,2*n_gauss+2)+sca*wk(ik)*exp((et(jj,ik)-e_fermi)*times(it))
               endif
            enddo

            !write(stdout,*) 'ex(', ii ,',1) = ', ex(ii,1)
            !call flush_unit(stdout)

            if(l_gauss .and. .not.l_truncated_coulomb) then
               if(ii==jj) then
                  do it=1,n_gauss+1
                     if(it/=n_gauss+1) then
                        ex(ii,it)=ex(ii,it)-wk(ik)*exp((et(jj,ik)-e_fermi)*times(it))*f_gauss(ik)
                     else
                        ex(ii,2*n_gauss+2)=ex(ii,2*n_gauss+2)-wk(ik)*exp((et(jj,ik)-e_fermi)*times(it))*f_gauss(ik)
                     endif
                  enddo
!add also the integral part divided by the number of k points
                  do it=1,n_gauss+1
                     if(it/=n_gauss+1) then
                        ex(ii,it)=ex(ii,it)+wk(ik)*exp((et(jj,ik)-e_fermi)*times(it))*f_gauss_int
                     else
                        ex(ii,2*n_gauss+2)=ex(ii,2*n_gauss+2)+wk(ik)*exp((et(jj,ik)-e_fermi)*times(it))*f_gauss_int
                     endif
                  enddo

               endif
            endif
            call stop_clock('k other')
         enddo
         !loop on v states
         do jj=1,nbnd_v


            !calcolate charge

            prods(:)=conjg(wfcs_s(:,ii))*wfcs_vc(:,jj)



            IF ( doublegrid ) THEN
               CALL cinterpolate( prod, prods, 1 )
            ELSE
               prod(:)=prods(:)
            ENDIF
!ULTRASOFT STUFF

            if(okvan) call adduspos_r(prod,bec_s(:,ii),bec_vc(:,jj))

            CALL fwfft ('Dense', prod, dfftp)


!calculate energy
            sca=0.d0
            do ig = 1, ngm_max
               sca=sca+conjg(prod(nl(ig)))*prod(nl(ig))*fac(ig)
            enddo
!            write(stdout,*) 'SCA', ii,jj,sca!ATTENZIONE
            call mp_sum(sca)
            !loop on positive times
            do it=n_gauss+1, 2*n_gauss+1
               ex(ii,it)=ex(ii,it)+sca*wk(ik)*(-1.d0)*exp((et(jj,ik)-e_fermi)*times(it))
            enddo

            if(l_gauss .and. .not.l_truncated_coulomb) then
               if(ii==jj) then
                  do it=n_gauss+1, 2*n_gauss+1
                     ex(ii,it)=ex(ii,it)-(-1.d0)*wk(ik)*exp((et(jj,ik)-e_fermi)*times(it))*f_gauss(ik)
                  enddo
!add also the integral part divided by the number of k points
                  do it=n_gauss+1, 2*n_gauss+1
                     ex(ii,it)=ex(ii,it)+(-1.d0)*wk(ik)*exp((et(jj,ik)-e_fermi)*times(it))*f_gauss_int
                  enddo
               endif
            endif


         enddo
         call print_clock('k prod')
         call print_clock('k inte')
         call print_clock('k add')
         call print_clock('k fft')
         call print_clock('k sum')
         call print_clock('k other')

      enddo

   enddo

   ex(:,:)=ex(:,:)/2.d0!for the spin

!write perturbative HF energies on screen
   do ii=1,nbnd_s
      write(stdout,*) 'X energy', ii, ex(ii,n_gauss+1)!*rytoev
   enddo
   do it=1,2*n_gauss+2
      write(stdout,*) it,ex(nbnd_v+1,it),ex(nbnd_v+2,it)
   enddo


!write on file exchange parte

   if (ionode) then
      iun = find_free_unit()
      open(unit=iun,file=trim(prefix)//'.exchange',status='unknown',form='unformatted')
      write(iun) nbnd_s
      write(iun) ex(1:nbnd_s,n_gauss+1)
      close(iun)
   endif
!write on file time dependent part

   if (ionode) then
      iun = find_free_unit()
      open(unit=iun,file=trim(prefix)//'.gv_time',status='unknown',form='unformatted')
      write(iun) nbnd_s
      write(iun) n_gauss
      write(iun) tau_gauss
      do it=1,2*n_gauss+2
         write(iun) ex(1:nbnd_s,it)
      enddo
      close(iun)
   endif


   call print_clock('k prod')
   call print_clock('k inte')
   call print_clock('k add')
   call print_clock('k fft')
   call print_clock('k sum')
   call print_clock('k other')

   deallocate(wfcs_s, wfcs_vc)
   deallocate(bec_vc, bec_s)
   deallocate(fac)
   deallocate(prods,prod)
   deallocate(ex)
   deallocate(times,x,w)
   deallocate(f_gauss)

! #endif

 end subroutine dft_exchange_k


subroutine dft_exchange(nbnd_v,nbnd_s,n_set)
!this subroutine calculates the exchange
!energi term for every state and writes on disk

! #ifdef __GWW

  USE io_global,            ONLY : stdout, ionode, ionode_id
  USE io_files,             ONLY : find_free_unit, prefix, iunwfc, nwordwfc, iunigk
  USE mp_global,            ONLY : nproc_pool, me_pool
  USE kinds,                ONLY : DP
  USE basis
  USE klist
  USE constants,            ONLY : e2, pi, tpi, fpi, rytoev
  USE wvfct,                ONLY : igk, g2kin, npwx, npw, nbnd, nbndx, ecutwfc
  USE control_flags,        ONLY: gamma_only
  USE cell_base,            ONLY: at, alat, tpiba, omega, tpiba2
  USE wannier_gw
  USE fft_base,             ONLY : dffts, dfftp
  USE fft_interfaces,       ONLY : fwfft, invfft
  USE gvect
  USE gvecs,              ONLY : nls, nlsm, doublegrid
  USE uspp
  USE wavefunctions_module, ONLY : psic, evc
  USE realus,               ONLY : adduspos_gamma_r
  USE cell_base,            ONLY : at, bg, omega
  USE mp,                   ONLY : mp_sum, mp_bcast
  USE exx,                  ONLY : exx_divergence, exx_grid_init, yukawa
  USE becmod,               ONLY : calbec

  implicit none

   INTEGER, INTENT(in)  :: nbnd_v !number of valence states
   INTEGER, INTENT(in)  :: nbnd_s !number of states considered
   INTEGER, INTENT(in)  :: n_set  !defines the number of states to be read from disk at the same time
   REAL(kind=DP), ALLOCATABLE :: fac(:)
   REAL(kind=DP), ALLOCATABLE :: e_x(:)
   REAL(kind=DP) :: qq_fact,exxdiv
   INTEGER :: ig,iiv,iv,jjs,js,hw
   REAL(kind=DP), ALLOCATABLE :: becpr(:,:)
   REAL(kind=DP), ALLOCATABLE :: tmpreal1(:), tmpreal_v(:,:),tmpreal_s(:,:)
   INTEGER :: igk0(npwx)
   REAL(kind=dp) :: g2kin_bp(npwx)
   INTEGER :: npw0
   INTEGER :: jmin,jmax
   COMPLEX(kind=DP), ALLOCATABLE :: prod_g(:),prod_c(:)
   REAL(kind=DP), ALLOCATABLE :: prod_r(:)
   REAL(kind=DP) :: exc
   INTEGER :: iun

   allocate(fac(ngm))
!sets factors terms
!sets factors terms
!this has already  been called   call exx_grid_init()
   exxdiv=exx_divergence()
   do ig=1,ngm
      qq_fact = g(1,ig)**2.d0 + g(2,ig)**2.d0 + g(3,ig)**2.d0

      if(.not.l_truncated_coulomb) then

         if (qq_fact > 1.d-8) then
            fac(ig)=e2*fpi/(tpiba2*qq_fact + yukawa )
         else
            fac(ig)= - exxdiv ! & ! or rather something else (see F.Gygi)
            !                         - e2*fpi   ! THIS ONLY APPLYS TO HYDROGEN
            if (yukawa .gt. 1.d-8 ) then
               fac(ig) = fac(ig) + e2*fpi/(tpiba2*qq_fact + yukawa )
            end if
         end if
      else
         if (qq_fact > 1.d-8) then
            fac(ig)=(e2*fpi/(tpiba2*qq_fact))*(1.d0-dcos(dsqrt(qq_fact)*truncation_radius*tpiba))
         else
            fac(ig)=e2*fpi*(truncation_radius**2.d0/2.d0)

         endif
      endif

   end do
   fac(:)=fac(:)/omega

!   write(stdout,*) "truncation_radius", truncation_radius
!   write(stdout,*) "------------------------"
!   write(stdout,*) "fac()", fac, omega
!   write(stdout,*) "------------------------"

   allocate(e_x(nbnd_s))
   e_x(:)=0.d0
   CALL gk_sort(xk(1,1),ngm,g,ecutwfc/tpiba2, &
              &    npw0,igk0,g2kin_bp)

!    write(stdout,*) "------------------------"
!    write(stdout,*)"WHO KNOWS npw0, g2kin_bp, igk0", npw0, g2kin_bp, igk0
!    write(stdout,*) "------------------------"

!   write(stdout,*) "------------------------"
!   write(stdout,*) "g()", g
!   write(stdout,*) "------------------------"

   if(okvan) allocate(becpr(nkb,nbnd_s))
   if  ( nkb > 0 .and. okvan) then
      CALL init_us_2( npw, igk, xk(1,1), vkb )
      CALL calbec(npw, vkb, evc, becpr, nbnd_s)
   endif


   allocate(tmpreal1(dffts%nnr))
   allocate(tmpreal_v(dfftp%nnr,n_set))
   allocate(tmpreal_s(dfftp%nnr,n_set))
   allocate(prod_g(ngm))
   allocate(prod_c(dfftp%nnr))
   allocate(prod_r(dfftp%nnr))

!external loop on valence state
   do iiv=1,ceiling(real(nbnd_v)/real(n_set))
   !read states and do fourier transform
      do hw=(iiv-1)*n_set+1,min(iiv*n_set,nbnd_v),2
         psic(:)=(0.d0,0.d0)
         psic(:)=(0.d0,0.d0)
         IF ( hw < min(iiv*n_set,nbnd_v)) then
            psic(nls(igk(1:npw0)))  = evc(1:npw0,hw) + &
                 ( 0.D0, 1.D0 ) * evc(1:npw0,hw+1)
            psic(nlsm(igk(1:npw0))) = CONJG( evc(1:npw,hw) - &
                 ( 0.D0, 1.D0 ) * evc(1:npw0,hw+1) )
         ELSE
            psic(nls(igk(1:npw0)))  = evc(1:npw0,hw)
            psic(nlsm(igk(1:npw0))) = CONJG( evc(1:npw0,hw) )
         END IF

         CALL invfft ('Wave', psic, dffts)
         tmpreal1(1:dffts%nnr)=dble(psic(1:dffts%nnr))
         if(doublegrid) then
            call interpolate(tmpreal_v(:,hw-(iiv-1)*n_set),tmpreal1,1)
         else
          tmpreal_v(:,hw-(iiv-1)*n_set)=tmpreal1(:)
       endif
       if ( hw < min(iiv*n_set,nbnd_v)) then
          tmpreal1(1:dffts%nnr)=aimag(psic(1:dffts%nnr))
          if(doublegrid) then
             call interpolate(tmpreal_v(:,hw-(iiv-1)*n_set+1),tmpreal1,1)
          else
             tmpreal_v(:,hw-(iiv-1)*n_set+1)=tmpreal1(:)
          endif
       endif

    enddo

    do jjs=iiv,ceiling(real(nbnd_s)/real(n_set))
   !external loop on states
      !read states and do fourier transform
       do hw=(jjs-1)*n_set+1,min(jjs*n_set,nbnd_s),2
         psic(:)=(0.d0,0.d0)
         psic(:)=(0.d0,0.d0)
         IF ( hw < min(jjs*n_set,nbnd_s)) then
            psic(nls(igk(1:npw0)))  = evc(1:npw0,hw) + &
                 ( 0.D0, 1.D0 ) * evc(1:npw0,hw+1)
            psic(nlsm(igk(1:npw0))) = CONJG( evc(1:npw,hw) - &
                 ( 0.D0, 1.D0 ) * evc(1:npw0,hw+1) )
         ELSE
            psic(nls(igk(1:npw0)))  = evc(1:npw0,hw)
            psic(nlsm(igk(1:npw0))) = CONJG( evc(1:npw0,hw) )
         END IF

         CALL invfft ('Wave', psic, dffts)
         tmpreal1(1:dffts%nnr)=dble(psic(1:dffts%nnr))
         if(doublegrid) then
            call interpolate(tmpreal_s(:,hw-(jjs-1)*n_set),tmpreal1,1)
         else
          tmpreal_s(:,hw-(jjs-1)*n_set)=tmpreal1(:)
       endif
       if ( hw < min(jjs*n_set,nbnd_s)) then
          tmpreal1(1:dffts%nnr)=aimag(psic(1:dffts%nnr))
          if(doublegrid) then
             call interpolate(tmpreal_s(:,hw-(jjs-1)*n_set+1),tmpreal1,1)
          else
             tmpreal_s(:,hw-(jjs-1)*n_set+1)=tmpreal1(:)
          endif
       endif

    enddo

      !internal loop on valence states
         do iv=(iiv-1)*n_set+1,min(iiv*n_set,nbnd_v)
            !internal loop on states
            if(iiv==jjs) then
               jmin=iv
!                jmin=(jjs-1)*n_set+1
            else
               jmin=(jjs-1)*n_set+1
            endif
            jmax=min(jjs*n_set,nbnd_s)

            do js=jmin,jmax
               !do product in real speace
               prod_r(:)=tmpreal_v(:,iv-(iiv-1)*n_set)*tmpreal_s(:,js-(jjs-1)*n_set)
               if(okvan) call adduspos_gamma_r &
               (iv,js, prod_r(:),1,becpr(:,iv),becpr(:,js))

               prod_c(:)=dcmplx(prod_r(:),0.d0)
               CALL fwfft ('Dense', prod_c, dfftp)
               !go to g_space
               prod_g(1:ngm)=prod_c(nl(1:ngm))
               !calculated exchange
               exc=0.d0
               do ig=1,ngm
                  exc=exc+2.d0*dble(conjg(prod_g(ig))*prod_g(ig))*fac(ig)
               enddo
               if(gstart==2) exc=exc-dble(prod_g(1))*dble(prod_g(1))*fac(1)
               call mp_sum(exc)
               exc=-exc
               e_x(js)=e_x(js)+exc
               if(js<=nbnd_v .and. js /= iv)  then
                  e_x(iv)=e_x(iv)+exc
               endif
            enddo
         enddo
      enddo
   enddo

   do iv=1,nbnd_s
      write(stdout,*) 'Exchange energy', iv, e_x(iv)
   enddo
!write on file

   if (ionode) then
      iun = find_free_unit()
      open(unit=iun,file=trim(prefix)//'.exchange',status='unknown',form='unformatted')
      write(iun) nbnd_s
      write(iun) e_x(1:nbnd_s)
      close(iun)
   endif




   deallocate(tmpreal1,tmpreal_s,tmpreal_v)

   deallocate(fac)
   deallocate(e_x)
   deallocate(prod_c,prod_g)
   deallocate(prod_r)
   if(okvan) deallocate(becpr)

! #endif
 end subroutine dft_exchange


!----------------------------------------------------------------------
subroutine addus_charge(r_ij,becp_iw,becp_jw)
  !----------------------------------------------------------------------
  !
  !  This routine adds to the charge density the part which is due to
  !  the US augmentation.
  !
! #ifdef __GWW
  !
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat, ntyp => nsp, ityp
  USE fft_base,             ONLY : dfftp
  USE fft_interfaces,       ONLY : fwfft, invfft
  USE gvect,                ONLY : ngm, nl, nlm, gg, g, &
                                   eigts1, eigts2, eigts3, mill
  USE lsda_mod,             ONLY : nspin
  USE scf,                  ONLY : rho
  USE uspp,                 ONLY : okvan, nkb
  USE uspp_param,           ONLY : upf, lmaxq, nh
  USE control_flags,        ONLY: gamma_only
  USE wavefunctions_module, ONLY : psic
  !
  implicit none

  COMPLEX(kind=DP), INTENT(inout) :: r_ij(dfftp%nnr)!where to add the us term
  COMPLEX(kind=DP), INTENT(in) ::  becp_iw( nkb)!overlap of wfcs with us  projectors
  COMPLEX(kind=DP), INTENT(in) ::  becp_jw( nkb)!overlap of wfcs with us  projectors



  !
  !     here the local variables
  !

  integer :: ig, na, nt, ih, jh, is
  ! counters

  real(DP), allocatable :: qmod (:), ylmk0 (:,:)
  ! the modulus of G
  ! the spherical harmonics

  complex(DP) :: skk
  complex(DP), allocatable ::  aux (:,:), qgm(:)
  ! work space for rho(G,nspin)
  ! Fourier transform of q
  INTEGER, ALLOCATABLE :: ind_cor(:,:,:)
  INTEGER :: ijkb0, ikb,np

  if (.not.okvan) return


  allocate (aux ( ngm, nspin))
  allocate (qmod( ngm))
  allocate (qgm( ngm))
  allocate (ylmk0( ngm, lmaxq * lmaxq))

  aux (:,:) = (0.d0, 0.d0)
  call ylmr2 (lmaxq * lmaxq, ngm, g, gg, ylmk0)
  do ig = 1, ngm
     qmod (ig) = sqrt (gg (ig) )
  enddo

!found index correspondence


  allocate(ind_cor(ntyp,nat,maxval(nh(1:ntyp))))

  ijkb0 = 0
  do  np = 1, ntyp
     if ( upf(np)%tvanp ) then
        do  na = 1, nat
           if ( ityp(na) == np ) then
              do ih = 1, nh(np)
                ikb = ijkb0 + ih
                ind_cor(np,na,ih)=ikb
             enddo
             ijkb0=ijkb0+nh(np)
          endif
       enddo
    else
       do na=1,nat
          if(ityp(na) == np) ijkb0=ijkb0+nh(np)
       enddo
    endif
 enddo








  do nt = 1, ntyp
     if ( upf(nt)%tvanp ) then
        do ih = 1, nh (nt)
           do jh = 1, nh (nt)
              call qvan2 (ngm, ih, jh, nt, qmod, qgm, ylmk0)
              do na = 1, nat
                 if (ityp (na) .eq.nt) then
                    !
                    !  Multiply becsum and qg with the correct structure factor
                    !
                    do is = 1, nspin
                       do ig = 1, ngm
                          skk = eigts1 (mill (1,ig), na) * &
                                eigts2 (mill (2,ig), na) * &
                                eigts3 (mill (3,ig), na)
                          aux(ig,is)=aux(ig,is) + qgm(ig)*skk*&
                                  &conjg(becp_iw(ind_cor(nt,na,ih)))*becp_jw(ind_cor(nt,na,jh))
                       enddo
                    enddo
                 endif
              enddo
           enddo
        enddo
     endif
  enddo

  deallocate(ind_cor)
  !
  deallocate (ylmk0)
  deallocate (qgm)
  deallocate (qmod)
  !
  !     convert aux to real space and add to the charge density
  !
  do is = 1, nspin!SPIN TO BE IMPLEMENTED YET
     psic(:) = (0.d0, 0.d0)
     psic( nl(:) ) = aux(:,is)
     if (gamma_only) psic( nlm(:) ) = CONJG(aux(:,is))
     CALL invfft ('Dense', psic, dfftp)
     r_ij(:)=r_ij(:)+psic(:)
  enddo
  deallocate (aux)

! #endif
  return
end subroutine addus_charge

