! FOR GWW
!
! Author: P. Umari
! Modified by G. Stenuit
!
 subroutine create_vcprim(n_set, lzero, orthonorm,ecutoff,lcprim)
!this subroutine calculates the elemnents S_c',c,i=\int dr Psi_c'(r) Psi_c(r) \tilde{w^P_i}
!and writes on disk for each c'

! #ifdef __GWW
   USE io_global,            ONLY : stdout, ionode
   USE io_files,             ONLY : find_free_unit, prefix, diropn
   USE kinds,                ONLY : DP
   USE wannier_gw
   USE gvect
   USE constants,            ONLY : e2, pi, tpi, fpi
   USE cell_base,            ONLY : at, alat, tpiba, omega, tpiba2
   USE exx,                  ONLY : exx_divergence, yukawa
   USE wvfct,                ONLY : igk, g2kin, npwx, npw, nbnd, nbndx
   USE wavefunctions_module, ONLY : evc, psic
   USE fft_base,             ONLY : dffts, dfftp
   USE fft_interfaces,       ONLY : fwfft, invfft
   USE gvecs,              ONLY : nls, nlsm, doublegrid
   USE uspp,                 ONLY : nkb, vkb, okvan
   USE realus,               ONLY : adduspos_gamma_r
   USE mp,                   ONLY : mp_sum
   USE mp_global,            ONLY : mpime
   USE klist,                ONLY : xk
   USE memory_limit
   USE becmod,               ONLY : calbec

  implicit none

  INTEGER, INTENT(in)  :: n_set  !defines the number of product states to be read from disk at the same time
  LOGICAL, INTENT(in)  :: lzero !if true put the term G=0,G=0 of v to zero
  INTEGER, INTENT(in)  :: orthonorm!if == 1 opens orthonormalized products of wannier file, if==2 opens reduced orthonormalized products
  REAL(kind=DP), INTENT(in) :: ecutoff!cutoff in Rydberg for g sum
  LOGICAL,INTENT(in) :: lcprim!if true calculates only cprim terms

  REAL(kind=DP) :: exxdiv, qq
  INTEGER :: ig
  REAL(kind=DP), ALLOCATABLE :: fac(:)
  INTEGER :: ii, jj, iw, ir
  REAL(kind=DP), ALLOCATABLE :: smat(:,:)
  INTEGER :: iunsterms, iungprod
  CHARACTER(4) :: nfile
  REAL(kind=DP), ALLOCATABLE :: tmp_reals(:,:), tmp_prod(:),tmp_reals_jj(:)
  COMPLEX(kind=DP), ALLOCATABLE :: tmp_g(:), tmp_w(:),  tmp_wp(:,:)
  LOGICAL :: exst
  REAL(kind=DP) :: sca
  INTEGER :: n_buf,ib,i_first,i_last,idumm
  INTEGER :: ngm_max
  INTEGER :: nbnd_first
  REAL(kind=DP), ALLOCATABLE :: sca_vec(:)
  INTEGER :: n_set_new  !! NEW variables

  write(stdout,*) 'Routine create_cprim'
  call flush_unit(stdout)

!determine ngm_max
  ngm_max=0
  do ig=1,ngm
     if(gg(ig)*tpiba2 >= ecutoff) exit
     ngm_max=ngm_max+1
  enddo

  write(stdout,*) 'NGM MAX:', ngm_max, ngm
  !! NEW PART
  !!! evaluation of the new n_set: n_set_new to optimize the speed vs memry used:

  mem_used=0
  mem_used=mem_used+ngm_max*8
  mem_used=mem_used+dfftp%nnr*8
  mem_used=mem_used+numw_prod*nbnd_normal*8
  mem_used=mem_used+dffts%nnr*(cprim_last-cprim_first+1)*8
  mem_used=mem_used+dfftp%nnr*8
  mem_used=mem_used+ngm*16

!calculate becs

  if  ( nkb > 0 .and. okvan) then
     CALL init_us_2( npw, igk, xk(1,1), vkb )
     CALL calbec(npw, vkb, evc, becp_gw_c, nbnd)
  endif


  if(.not.lzero .and. .not.l_truncated_coulomb) then
     exxdiv=exx_divergence()
  else
     exxdiv = 0.d0
  endif

  write(stdout,*) 'ATTENZIONE1'
  call flush_unit(stdout)

  !calculate V(G)
  allocate(fac(ngm_max))
  allocate(tmp_prod(dfftp%nnr))
  !allocate(sca_vec(n_set_new))


  do ig=1,ngm_max
     qq = g(1,ig)**2.d0 + g(2,ig)**2.d0 + g(3,ig)**2.d0

     if(.not.l_truncated_coulomb) then

        if (qq > 1.d-8) then
           fac(ig)=e2*fpi/(tpiba2*qq + yukawa )
        else
           fac(ig)= - exxdiv
           if (yukawa .gt. 1.d-8 ) then
              fac(ig) = fac(ig) + e2*fpi/(tpiba2*qq + yukawa )
           end if
        end if
     else
        if (qq > 1.d-8) then
           fac(ig)=(e2*fpi/(tpiba2*qq))*(1.d0-dcos(dsqrt(qq)*truncation_radius*tpiba))
        else
           fac(ig)=e2*fpi*(truncation_radius**2.d0/2.d0)
        endif
     endif

  end do
    write(stdout,*) 'ATTENZIONE2'
    call flush_unit(stdout)

  fac(:)=fac(:)/omega
  if(lzero .and. gstart == 2) fac(1)=0.d0


  allocate( smat(numw_prod,nbnd_normal))
  allocate(tmp_reals(dffts%nnr, cprim_last-cprim_first+1),tmp_reals_jj(dfftp%nnr))
  allocate(tmp_g(ngm))

  write(stdout,*) 'ATTENZIONE3'
  write(stdout,*) 'lbound and ubound of psic: ', lbound(psic), ubound(psic)
  call flush_unit(stdout)

!put states on reals grid
  do ii = cprim_first, cprim_last, 2
     psic(:) = ( 0.D0, 0.D0 )
     IF ( ii < cprim_last) THEN
        psic(nls(igk(1:npw)))  = evc(1:npw,ii) + &
             ( 0.D0, 1.D0 ) * evc(1:npw,ii+1)
        psic(nlsm(igk(1:npw))) = CONJG( evc(1:npw,ii) - &
             ( 0.D0, 1.D0 ) * evc(1:npw,ii+1) )
     ELSE
        psic(nls(igk(1:npw)))  = evc(1:npw,ii)
        psic(nlsm(igk(1:npw))) = CONJG( evc(1:npw,ii) )
     END IF
     CALL invfft ('Wave', psic, dffts)
     tmp_reals(:, ii-cprim_first+1)= DBLE(psic(:))
     if(ii+1 <= cprim_last) then
        tmp_reals(:, ii+1-cprim_first+1)=dimag(psic(:))
     endif
  enddo

  write(stdout,*) 'ATTENZIONE4'
  call flush_unit(stdout)

!open product of wanniers filed
   iungprod = find_free_unit()
   if(orthonorm==0) then
      CALL diropn( iungprod, 'wiwjwfc', max_ngm*2, exst )
   else if(orthonorm==1) then
      CALL diropn( iungprod, 'wiwjwfc_on', max_ngm*2, exst )
   else
      CALL diropn( iungprod, 'wiwjwfc_red', max_ngm*2, exst )
   endif

!loop on w_p bufferized on n_set_new
  !ncpu=number_cpu()
  n_set_new=numw_prod
  call test_allocation_mp(ngm_max, n_set_new, 2, mem_used)


  write(stdout,*) '=> n_set_new=', n_set_new
  call flush_unit(stdout)

  n_buf=numw_prod/n_set_new
  if(n_buf*n_set_new < numw_prod) n_buf=n_buf+1

  write(stdout,*) '=> n_buf=', n_buf
  call flush_unit(stdout)
  allocate(tmp_wp(ngm_max,n_set_new))
  allocate(sca_vec(n_set_new))

  write(stdout,*) 'ATTENZIONE5'
  call flush_unit(stdout)


  do ib=1,n_buf

     i_first=(ib-1)*n_set_new+1
     i_last=min(ib*n_set_new,numw_prod)

!read in products of wanniers in G space and trasform to dense R space
     do iw=i_first,i_last
        call davcio(tmp_g,max_ngm*2,iungprod,iw,-1)
        tmp_wp(1:ngm_max,iw-i_first+1)=tmp_g(1:ngm_max)
     enddo

!!!!!!!!!
!loop on c'
     do ii=cprim_first,cprim_last
        write(stdout,*) 'State:', ii
        call flush_unit(stdout)
 !open file S_c'

!N in case re-read from disk


        if(ionode) then
           iunsterms =  find_free_unit()
           write(nfile,'(4i1)') &
                & ii/1000,mod(ii,1000)/100,mod(ii,100)/10,mod(ii,10)

           if(ib==1) then
              if(.not.lcprim) then
                 open( unit= iunsterms, file=trim(prefix)//'.vcprim.'//nfile, status='unknown',form='unformatted')
              else
                 open( unit= iunsterms, file=trim(prefix)//'.cprim.'//nfile, status='unknown',form='unformatted')
              endif
              smat(:,:)=0.d0
           else
              if(.not.lcprim) then
                 open( unit= iunsterms, file=trim(prefix)//'.vcprim.'//nfile, status='old',form='unformatted')
              else
                 open( unit= iunsterms, file=trim(prefix)//'.cprim.'//nfile, status='old',form='unformatted')
              endif
              read(iunsterms) idumm
              read(iunsterms) idumm
              read(iunsterms) idumm
              read(iunsterms) idumm
              do jj=1,nbnd_normal
                 read(iunsterms) smat(1:numw_prod, jj)
              enddo
              rewind(iunsterms)
           endif
        endif


 !loop on c
        if(.not.lcprim) then
           nbnd_first=1
        else
           nbnd_first=num_nbndv+1
        endif

        do jj=nbnd_first,nbnd_normal
           psic(:) = ( 0.D0, 0.D0 )
           psic(nls(igk(1:npw)))  = evc(1:npw,jj)
           psic(nlsm(igk(1:npw))) = CONJG( evc(1:npw,jj) )

           CALL invfft ('Wave', psic, dffts)
           tmp_reals_jj(:)= DBLE(psic(:))

           tmp_prod(1:dffts%nnr)=tmp_reals(:,ii-cprim_first+1)*tmp_reals_jj(:)
           if(doublegrid) then
              call interpolate(tmp_prod,tmp_prod,1)
           endif
           if(okvan) call adduspos_gamma_r(ii,jj,tmp_prod,1,becp_gw_c(:,ii),becp_gw_c(:,jj))
!trasform to g-space
           psic(1:dfftp%nnr)=cmplx(tmp_prod(1:dfftp%nnr),0.d0)
           CALL fwfft ('Dense', psic, dfftp)
           do ig=1,ngm_max
              tmp_g(ig)=psic(nl(ig))*fac(ig)
           enddo
!the following line for calling blas routine
           if(gstart==2) tmp_g(1) = 0.5d0*tmp_g(1)

           call dgemv('T',2*ngm_max,i_last-i_first+1,2.d0,tmp_wp,2*ngm_max,tmp_g,1,0.d0,sca_vec,1)

           call mp_sum(sca_vec(1: i_last-i_first+1))
           smat(i_first:i_last,jj)=sca_vec(1: i_last-i_first+1)
        enddo
  !write on file
        if(ionode) then
           write(iunsterms) ii
           write(iunsterms) num_nbndv
           write(iunsterms) nbnd_normal
           write(iunsterms) numw_prod
           do jj=nbnd_first,nbnd_normal
              write(iunsterms) smat(1:numw_prod, jj)
           enddo
           close(iunsterms)
        endif
     enddo

  enddo!on ib
  close(iungprod)
  deallocate(fac,smat, tmp_reals, tmp_prod,tmp_reals_jj)
  deallocate(tmp_g)
  deallocate(tmp_wp)
  deallocate(sca_vec)
  return
! #endif
end subroutine create_vcprim
 !
 !
subroutine create_vcw_overlap(n_set, orthonorm,ecutoff)
!this subroutine calculates the elemnents S_c',c,i=\int dr Psi_c'(r) Psi_c(r) \tilde{w^P_i}
!and writes on disk for each c'

! #ifdef __GWW
   USE io_global,            ONLY : stdout, ionode
   USE io_files,             ONLY : find_free_unit, prefix, diropn
   USE kinds,                ONLY : DP
   USE wannier_gw
   USE gvect
   USE constants,            ONLY : e2, pi, tpi, fpi
   USE cell_base,            ONLY : at, alat, tpiba, omega, tpiba2
   USE exx,                  ONLY : exx_divergence, yukawa
   USE wvfct,                ONLY : igk, g2kin, npwx, npw, nbnd, nbndx
   USE wavefunctions_module, ONLY : evc, psic
   USE fft_base,             ONLY : dffts, dfftp
   USE fft_interfaces,       ONLY : fwfft, invfft
   USE gvecs,              ONLY : nls, nlsm, doublegrid
   USE uspp,                 ONLY : nkb, vkb, okvan
   USE realus,               ONLY : adduspos_gamma_r
   USE mp,                   ONLY : mp_sum
   USE mp_global,            ONLY : mpime
   USE klist,                ONLY : xk
   USE memory_limit
   USE becmod,               ONLY : calbec

  implicit none

  INTEGER, INTENT(in)  :: n_set  !defines the number of states to be read from disk at the same time
  INTEGER, INTENT(in)  :: orthonorm!if true opens orthonormalized products of wannier file
  REAL(kind=DP), INTENT(in) :: ecutoff!cutoff in Rydberg for g sum

  INTEGER :: ig
  INTEGER :: ii, jj, iw, ir
  REAL(kind=DP), ALLOCATABLE :: smat(:,:)
  INTEGER :: iunsterms, iungprod
  CHARACTER(4) :: nfile
  REAL(kind=DP), ALLOCATABLE :: tmp_reals(:,:), tmp_prod(:),tmp_reals_jj(:)
  COMPLEX(kind=DP), ALLOCATABLE :: tmp_g(:), tmp_w(:), tmp_wp(:,:)
  LOGICAL :: exst
  REAL(kind=DP) :: sca
  INTEGER :: n_buf,ib,i_first,i_last,idumm
  INTEGER :: ngm_max
  REAL(kind=DP), ALLOCATABLE :: sca_vec(:)
  INTEGER :: cprim_last_eff
  INTEGER :: n_set_new  !! NEW variables

  write(stdout,*) 'Routine create_vcw_overlap'
  call flush_unit(stdout)


  cprim_last_eff=min(cprim_last,num_nbndv)

!determine ngm_max
  ngm_max=0
  do ig=1,ngm
     if(gg(ig)*tpiba2 >= ecutoff) exit
     ngm_max=ngm_max+1
  enddo
  write(stdout,*) 'NGM MAX:', ngm_max, ngm

  mem_used=0
  mem_used=mem_used+dfftp%nnr*8
  mem_used=mem_used+numw_prod*(nbnd_normal-num_nbndv)*8
  mem_used=mem_used+dffts%nnr*(cprim_last_eff-cprim_first+1)*8
  mem_used=mem_used+dffts%nnr*8
  mem_used=mem_used+ngm*16

!calculate becs

  if  ( nkb > 0 .and. okvan) then
     CALL init_us_2( npw, igk, xk(1,1), vkb )
     CALL calbec(npw, vkb, evc, becp_gw_c, nbnd)
  endif



  write(stdout,*) 'ATTENZIONE1'
  call flush_unit(stdout)

  !calculate V(G)
  allocate(tmp_prod(dfftp%nnr))

  allocate( smat(numw_prod,nbnd_normal-num_nbndv))
  allocate(tmp_reals(dffts%nnr, cprim_last_eff-cprim_first+1),tmp_reals_jj(dffts%nnr))
  allocate(tmp_g(ngm))
  ! allocate(sca_vec(n_set_new))
  write(stdout,*) 'ATTENZIONE3'
  call flush_unit(stdout)

!put states on reals grid
  do ii = cprim_first, cprim_last_eff, 2
     psic(:) = ( 0.D0, 0.D0 )
     IF ( ii < cprim_last_eff ) THEN
        psic(nls(igk(1:npw)))  = evc(1:npw,ii) + &
             ( 0.D0, 1.D0 ) * evc(1:npw,ii+1)
        psic(nlsm(igk(1:npw))) = CONJG( evc(1:npw,ii) - &
             ( 0.D0, 1.D0 ) * evc(1:npw,ii+1) )
     ELSE
        psic(nls(igk(1:npw)))  = evc(1:npw,ii)
        psic(nlsm(igk(1:npw))) = CONJG( evc(1:npw,ii) )
     END IF
     CALL invfft ('Wave', psic, dffts)
     tmp_reals(:, ii-cprim_first+1)= DBLE(psic(:))
     if(ii+1 <= cprim_last_eff) then
        tmp_reals(:, ii+1-cprim_first+1)=dimag(psic(:))
     endif
  enddo

  write(stdout,*) 'ATTENZIONE4'
  call flush_unit(stdout)

!open product of wanniers filed
   iungprod = find_free_unit()
   if(orthonorm == 0) then
      CALL diropn( iungprod, 'wiwjwfc', max_ngm*2, exst )
   else if(orthonorm==1) then
      CALL diropn( iungprod, 'wiwjwfc_on', max_ngm*2, exst )
   else
      CALL diropn( iungprod, 'wiwjwfc_red', max_ngm*2, exst )
   endif

!loop on w_p bufferized on n_set_new

  n_set_new=numw_prod
  call test_allocation_mp(ngm_max,n_set_new,2,mem_used)
  write(stdout,*) '=> n_set_new=', n_set_new
  call flush_unit(stdout)


  n_buf=numw_prod/n_set_new
  if(n_buf*n_set_new < numw_prod) n_buf=n_buf+1
  allocate(tmp_wp(ngm_max,n_set_new))
  allocate(sca_vec(n_set_new))

  write(stdout,*) 'After allocate(tmp_wp(ngm_max,n_set_new))'
  write(stdout,*) '=> n_buf=', n_buf
  call flush_unit(stdout)

  do ib=1,n_buf

     i_first=(ib-1)*n_set_new+1
     i_last=min(ib*n_set_new,numw_prod)

!read in products of wanniers in G space and trasform to dense R space
     do iw=i_first,i_last
        call davcio(tmp_g,max_ngm*2,iungprod,iw,-1)
        tmp_wp(1:ngm_max,iw-i_first+1)=tmp_g(1:ngm_max)
     enddo

!!!!!!!!!
!loop on c'
     do ii=cprim_first,min(cprim_last,num_nbndv)
        write(stdout,*) 'State:', ii
        call flush_unit(stdout)
 !open file S_c'


!N in case re-read from disk


        if(ionode) then
           iunsterms =  find_free_unit()
           write(nfile,'(4i1)') &
                & ii/1000,mod(ii,1000)/100,mod(ii,100)/10,mod(ii,10)

           if(ib==1) then
              open( unit= iunsterms, file=trim(prefix)//'.vcw_overlap.'//nfile, status='unknown',form='unformatted')
              smat(:,:)=0.d0
           else
              open( unit= iunsterms, file=trim(prefix)//'.vcw_overlap.'//nfile, status='old',form='unformatted')
              read(iunsterms) idumm
              read(iunsterms) idumm
              read(iunsterms) idumm
              read(iunsterms) idumm
              do jj=1,nbnd_normal-num_nbndv
                 read(iunsterms) smat(1:numw_prod, jj)
              enddo
              rewind(iunsterms)
           endif
        endif


 !loop on c
        do jj=num_nbndv+1,nbnd_normal
           psic(:) = ( 0.D0, 0.D0 )
           psic(nls(igk(1:npw)))  = evc(1:npw,jj)
           psic(nlsm(igk(1:npw))) = CONJG( evc(1:npw,jj) )

           CALL invfft ('Wave', psic, dffts)
           tmp_reals_jj(:)= DBLE(psic(:))

           tmp_prod(1:dffts%nnr)=tmp_reals(:,ii-cprim_first+1)*tmp_reals_jj(:)
           if(doublegrid) then
              call interpolate(tmp_prod,tmp_prod,1)
           endif
           if(okvan) call adduspos_gamma_r(ii,jj,tmp_prod,1,becp_gw_c(:,ii),becp_gw_c(:,jj))
!trasform to g-space
           psic(1:dfftp%nnr)=cmplx(tmp_prod(1:dfftp%nnr),0.d0)
           CALL fwfft ('Dense', psic, dfftp)
           do ig=1,max_ngm
               tmp_g(ig)=psic(nl(ig))
            enddo

!the following line for calling blas routine
            if(gstart==2) tmp_g(1) = 0.5d0*tmp_g(1)

            call dgemv('T',2*ngm_max,i_last-i_first+1,2.d0,tmp_wp,2*ngm_max,tmp_g,1,0.d0,sca_vec, 1)

            call mp_sum(sca_vec(1: i_last-i_first+1))
            smat(i_first:i_last,jj-num_nbndv)=sca_vec(1: i_last-i_first+1)

         enddo
  !write on file
        if(ionode) then
           write(iunsterms) ii
           write(iunsterms) num_nbndv
           write(iunsterms) nbnd_normal
           write(iunsterms) numw_prod
           do jj=num_nbndv+1,nbnd_normal
              write(iunsterms) smat(1:numw_prod, jj-num_nbndv)
           enddo
           close(iunsterms)
        endif
     enddo

  enddo!on ib
  close(iungprod)
  deallocate(smat, tmp_reals, tmp_prod,tmp_reals_jj)
  deallocate(tmp_g)
  deallocate(tmp_wp)
  deallocate(sca_vec)
  return
! #endif
end subroutine create_vcw_overlap
 !
 !
subroutine create_upper_states(n_set, lzero, orthonorm,ecutoff)
!this subroutine calculates the reduced upper states and their products
!with products of wannier

! #ifdef __GWW
   USE io_global,            ONLY : stdout, ionode, ionode_id
   USE io_files,             ONLY : find_free_unit, prefix, diropn
   USE kinds,                ONLY : DP
   USE wannier_gw
   USE gvect
   USE constants,            ONLY : e2, pi, tpi, fpi
   USE cell_base,            ONLY : at, alat, tpiba, omega, tpiba2
   USE exx,                  ONLY : exx_divergence, yukawa
   USE wvfct,                ONLY : igk, g2kin, npwx, npw, nbnd, nbndx, et
   USE wavefunctions_module, ONLY : evc, psic
   USE fft_base,             ONLY : dffts, dfftp
   USE fft_interfaces,       ONLY : fwfft, invfft
   USE gvecs,              ONLY : nls, nlsm, doublegrid
   USE uspp,                 ONLY : nkb, vkb, okvan
   USE realus,               ONLY : adduspos_gamma_r
   USE mp,                   ONLY : mp_sum, mp_bcast
   USE mp_global,            ONLY : mpime
   USE klist,                ONLY : xk
   USE becmod,               ONLY : calbec

  implicit none

  INTEGER, INTENT(in)  :: n_set  !defines the number of product states to be read from disk at the same time
  LOGICAL, INTENT(in)  :: lzero !if true put the term G=0,G=0 of v to zero
  INTEGER, INTENT(in)  :: orthonorm!if == 1 opens orthonormalized products of wannier file,
                                    !if==2 opens reduced orthonormalized products
  REAL(kind=DP), INTENT(in) :: ecutoff!cutoff in Rydberg for g sum


  INTEGER :: ngm_max, ig, ii, jj, iw, jw, id, i, j, m
  REAL(kind=DP), ALLOCATABLE :: fac(:), tmp_reals(:,:), tmp_reals_up(:,:), tmp_prod(:)
  REAL(kind=DP) :: qq, exxdiv, sca
  COMPLEX(kind=DP), ALLOCATABLE :: tmp_wp(:,:)
  INTEGER :: iungprod
  LOGICAL :: exst
  INTEGER :: ndelta, id_first, id_last, id_states, id_upper, n_upper_tot
  REAL(kind=DP), ALLOCATABLE :: o_mat(:,:), o_eigen(:), o_vector(:,:)
  COMPLEX(kind=DP), ALLOCATABLE :: tmp_g(:,:),  tmp_sigma(:,:), tmp_buf(:)
  INTEGER, ALLOCATABLE :: iwork(:),ifail(:)
  INTEGER :: lwork,info
  REAL(kind=DP), ALLOCATABLE :: work(:)
  REAL(kind=DP), ALLOCATABLE :: smat(:,:),vmat(:,:)
  CHARACTER(4) :: nfile
  INTEGER :: iunsterms
  REAL(kind=DP), ALLOCATABLE :: ene_reduced(:)

  write(stdout,*) 'Routine create_upper_states'
  call flush_unit(stdout)

  if(nbnd_normal >= nbnd) return

!determine ngm_max
  ngm_max=0
  do ig=1,ngm
     if(gg(ig)*tpiba2 >= ecutoff) exit
     ngm_max=ngm_max+1
  enddo

  write(stdout,*) 'NGM MAX:', ngm_max, ngm

  !calculate becs
  if  ( nkb > 0 .and. okvan) then
     CALL init_us_2( npw, igk, xk(1,1), vkb )
     CALL calbec(npw, vkb, evc, becp_gw_c, nbnd)
  endif

 if(.not.lzero .and. .not.l_truncated_coulomb) then
     exxdiv=exx_divergence()
  else
     exxdiv = 0.d0
  endif

  !calculate V(G)
  allocate(fac(ngm_max))
  do ig=1,ngm_max
     qq = g(1,ig)**2.d0 + g(2,ig)**2.d0 + g(3,ig)**2.d0

     if(.not.l_truncated_coulomb) then

        if (qq > 1.d-8) then
           fac(ig)=e2*fpi/(tpiba2*qq + yukawa )
        else
           fac(ig)= - exxdiv
           if (yukawa .gt. 1.d-8 ) then
              fac(ig) = fac(ig) + e2*fpi/(tpiba2*qq + yukawa )
           end if
        end if
     else
        if (qq > 1.d-8) then
           fac(ig)=(e2*fpi/(tpiba2*qq))*(1.d0-dcos(dsqrt(qq)*truncation_radius*tpiba))
        else
           fac(ig)=e2*fpi*(truncation_radius**2.d0/2.d0)
        endif
     endif

  end do
  fac(:)=fac(:)/omega
  if(lzero .and. gstart == 2) fac(1)=0.d0

  allocate(tmp_reals(dffts%nnr, cprim_last-cprim_first+1))

!put states on reals grid
  do ii = cprim_first, cprim_last, 2
     psic(:) = ( 0.D0, 0.D0 )
     IF ( ii < cprim_last) THEN
        psic(nls(igk(1:npw)))  = evc(1:npw,ii) + &
             ( 0.D0, 1.D0 ) * evc(1:npw,ii+1)
        psic(nlsm(igk(1:npw))) = CONJG( evc(1:npw,ii) - &
             ( 0.D0, 1.D0 ) * evc(1:npw,ii+1) )
     ELSE
        psic(nls(igk(1:npw)))  = evc(1:npw,ii)
        psic(nlsm(igk(1:npw))) = CONJG( evc(1:npw,ii) )
     END IF
     CALL invfft ('Wave', psic, dffts)
     tmp_reals(:, ii-cprim_first+1)= DBLE(psic(:))
     if(ii+1 <= cprim_last) then
        tmp_reals(:, ii+1-cprim_first+1)=dimag(psic(:))
     endif
  enddo


!if n_set >= numw_prod read in and store products of wanniers
!open product of wanniers filed
  if(n_set >= numw_prod) then
     iungprod = find_free_unit()
     if(orthonorm == 0) then
        CALL diropn( iungprod, 'wiwjwfc', max_ngm*2, exst )
     else if(orthonorm==1) then
        CALL diropn( iungprod, 'wiwjwfc_on', max_ngm*2, exst )
     else
        CALL diropn( iungprod, 'wiwjwfc_red', max_ngm*2, exst )
     endif

     allocate(tmp_wp(ngm_max,numw_prod))

!read in products of wanniers in G space
     allocate(tmp_buf(ngm))
     do iw=1,numw_prod
        call davcio(tmp_buf,max_ngm*2,iungprod,iw,-1)
        tmp_wp(1:ngm_max,iw)= tmp_buf(1:ngm_max)
     enddo
     deallocate(tmp_buf)
     close(iungprod)
  else
     allocate(tmp_wp(ngm,1))
  endif

!states are already in evc



!loop on delta
  ndelta=(nbnd-nbnd_normal)/num_nbnd_delta
  if(ndelta*num_nbnd_delta < (nbnd-nbnd_normal)) ndelta=ndelta+1
  allocate(tmp_reals_up(dffts%nnr, num_nbnd_delta))
  allocate( tmp_g(ngm_max, num_nbnd_delta))
  allocate(tmp_prod(dfftp%nnr))
  n_upper_tot=(ndelta-1)*num_nbnd_upper+min(nbnd-nbnd_normal-(ndelta-1)*num_nbnd_delta,num_nbnd_upper)
  allocate(ene_reduced(n_upper_tot))
  do id=1,ndelta
     write(stdout,*) 'ID:',id
     call flush_unit(stdout)
     id_first=nbnd_normal+(id-1)*num_nbnd_delta+1
     id_last=min(id_first+num_nbnd_delta-1,nbnd)
     id_states=id_last-id_first+1
     id_upper=min(id_states,num_nbnd_upper)
     write(stdout,*) 'ID:',id,id_first,id_last,id_states,id_upper
 !determined reduced energies
     sca=0.d0
     do ii=id_first,id_last
        sca=sca+et(ii,1)
     enddo
     sca=sca/dble(id_states)
     ene_reduced(1+(id-1)*num_nbnd_upper:1+(id-1)*num_nbnd_upper+id_upper-1)=sca

 !FFT states to R
     do ii = id_first, id_last, 2
        psic(:) = ( 0.D0, 0.D0 )
        IF ( ii < id_last) THEN
           psic(nls(igk(1:npw)))  = evc(1:npw,ii) + &
                ( 0.D0, 1.D0 ) * evc(1:npw,ii+1)
           psic(nlsm(igk(1:npw))) = CONJG( evc(1:npw,ii) - &
                ( 0.D0, 1.D0 ) * evc(1:npw,ii+1) )
        ELSE
           psic(nls(igk(1:npw)))  = evc(1:npw,ii)
           psic(nlsm(igk(1:npw))) = CONJG( evc(1:npw,ii) )
        END IF
        CALL invfft ('Wave', psic, dffts)
        tmp_reals_up(:, ii-id_first+1)= DBLE(psic(:))
        if(ii+1 <= id_last) then
           tmp_reals_up(:, ii+1-id_first+1)=dimag(psic(:))
        endif
     enddo
     allocate(o_mat(id_states,id_states))
     allocate(o_eigen(id_states), o_vector(id_states,id_upper))
    !loop on i
     do ii = cprim_first, cprim_last

       !form products
        do jj=1,id_states
           tmp_prod(1:dffts%nnr)=tmp_reals(:,ii-cprim_first+1)*tmp_reals_up(:,jj)
           if(doublegrid) then
              call interpolate(tmp_prod,tmp_prod,1)
           endif
           if(okvan) call adduspos_gamma_r &
                &(ii,jj,tmp_prod,1,becp_gw_c(:,ii),becp_gw_c(:,jj+id_first-1))
           psic(1:dfftp%nnr)=cmplx(tmp_prod(1:dfftp%nnr),0.d0)
           CALL fwfft ('Dense', psic, dfftp)
           do ig=1,ngm_max
              tmp_g(ig,jj)=psic(nl(ig))
           enddo
        enddo

       !calculate overlap matrix

    !    call zgemm('C','N',id_states,id_states,ngm_max,(1.d0,0.d0),tmp_g,&
    !            &ngm_max,tmp_g,ngm_max,(0.d0,0.d0),z_mat,id_states)
        call dgemm('T','N',id_states,id_states,2*ngm_max,2.d0,tmp_g,&
              &2*ngm_max,tmp_g,2*ngm_max,0.d0,o_mat,id_states)

        if(gstart==2) then
           do i=1,id_states
              do j=1,id_states
                 o_mat(j,i)=o_mat(j,i)-dble(conjg(tmp_g(1,j))*tmp_g(1,i))
              enddo
           enddo
        endif

        call mp_sum(o_mat(:,:))


        !calculate num_nbnd_upper upper eigenvalues/eigenvectors of o_mat
        if(mpime==0) then
           allocate(iwork(5*id_states), ifail(id_states))
           allocate(work(1))
           call dsyevx('V','I','U',id_states, o_mat, id_states,0.d0,0.d0,&
                &id_states-id_upper+1,id_states,0.d0,m,o_eigen,o_vector,id_states,&
                &work,-1,iwork,ifail,info)
           lwork=work(1)
           deallocate(work)
           allocate(work(lwork))
           call dsyevx('V','I','U',id_states, o_mat, id_states,0.d0,0.d0,&
                &id_states-id_upper+1,id_states,0.d0,m,o_eigen,o_vector,id_states,&
                &work,lwork,iwork,ifail,info)
           if(info/=0) then
              write(stdout,*) 'Routine create_upper_states: dsyevx failed'
              call flush_unit(stdout)
              stop
           endif
           deallocate(iwork,ifail,work)
        endif
        call mp_bcast(o_vector(:,:), ionode_id)
        call mp_bcast(o_eigen(:), ionode_id)
       !calculate reduced basis
        allocate( tmp_sigma(ngm_max,id_upper))

       !calculate reduced states
     !   call zgemm('N','N',ngm_max,id_upper,id_states,(1.d0,0.d0),tmp_g,ngm_max, z_vector,&
     !        &id_states, (0.d0,0.d0),tmp_sigma,ngm_max)
        call dgemm('N','N',2*ngm_max,id_upper,id_states,1.d0,tmp_g,2*ngm_max, o_vector,&
             &id_states, 0.d0,tmp_sigma,2*ngm_max)



       !calculate products
        allocate(vmat(numw_prod,id_upper),smat(numw_prod, id_upper))


        if(n_set >= numw_prod) then
        !overlap products
           if(ii <= num_nbndv) then

             ! call zgemm('C','N',numw_prod,id_upper,ngm_max,(1.d0,0.d0),tmp_wp,ngm_max, tmp_sigma,ngm_max, &
             !      &(0.d0,0.d0),z_smat,numw_prod)

              !smat(:,:)=2.d0*dble(z_smat(:,:))
               call dgemm('T','N',numw_prod,id_upper,2*ngm_max,2.d0,tmp_wp,2*ngm_max, tmp_sigma,2*ngm_max, &
                    &0.d0,smat,numw_prod)
              if(gstart==2) then
                 do iw=1,numw_prod
                    do jw=1,id_upper
                       smat(iw,jw)=smat(iw,jw)-dble(conjg(tmp_wp(1,iw))*tmp_sigma(1,jw))
                    enddo
                 enddo
              endif

           endif
           do iw=1,id_upper
              tmp_sigma(1:ngm_max,iw)=tmp_sigma(1:ngm_max,iw)*fac(1:ngm_max)
           enddo
!           call zgemm('C','N',numw_prod,id_upper,ngm_max,(1.d0,0.d0),tmp_wp,ngm_max, tmp_sigma,ngm_max, &
!                &(0.d0,0.d0),z_vmat,numw_prod)
!          vmat(:,:)=2.d0*dble(z_vmat(:,:))

           call dgemm('T','N',numw_prod,id_upper,2*ngm_max,2.d0,tmp_wp,2*ngm_max, tmp_sigma,2*ngm_max, &
                &0.d0,vmat,numw_prod)

           if(gstart==2) then
              do iw=1,numw_prod
                 do jw=1,id_upper
                    vmat(iw,jw)=vmat(iw,jw)-dble(conjg(tmp_wp(1,iw))*tmp_sigma(1,jw))
                 enddo
              enddo
           endif
        else
           iungprod = find_free_unit()
           if(orthonorm == 0) then
              CALL diropn( iungprod, 'wiwjwfc', max_ngm*2, exst )
           else if(orthonorm==1) then
              CALL diropn( iungprod, 'wiwjwfc_on', max_ngm*2, exst )
           else
              CALL diropn( iungprod, 'wiwjwfc_red', max_ngm*2, exst )
           endif
           do iw=1,numw_prod
              call davcio(tmp_wp(:,1),max_ngm*2,iungprod,iw,-1)
              do jw=1,id_upper
                 vmat(iw,jw)=0.d0
                 smat(iw,jw)=0.d0
                 do ig=1,ngm_max
                    vmat(iw,jw)=vmat(iw,jw)+2.d0*fac(ig)*dble(conjg(tmp_wp(ig,1))*tmp_sigma(ig,jw))
                 enddo
                 if(ii<=num_nbndv) then
                    do ig=1,ngm_max
                       smat(iw,jw)=smat(iw,jw)+2.d0*dble(conjg(tmp_wp(ig,1))*tmp_sigma(ig,jw))
                    enddo
                 endif
                 if(gstart==2) then
                    vmat(iw,jw)=vmat(iw,jw)-fac(1)*dble(conjg(tmp_wp(1,1))*tmp_sigma(1,jw))
                    if(ii<=num_nbndv) smat(iw,jw)=smat(iw,jw)-dble(conjg(tmp_wp(1,1))*tmp_sigma(1,jw))
                 endif
              enddo
           enddo
           close(iungprod)
        endif
        call mp_sum(vmat(:,:))

        if(ii<=num_nbndv) call mp_sum(smat(:,:))
 !write products on disk
        if(ionode) then
           iunsterms =  find_free_unit()
           write(nfile,'(4i1)') &
             & ii/1000,mod(ii,1000)/100,mod(ii,100)/10,mod(ii,10)
           !if there is the first block delta we start from the beginning
           if(ii<=num_nbndv) then
              if(id==1) then
                 open( unit= iunsterms, file=trim(prefix)//'.vcw_up_overlap.'//nfile, status='unknown',form='unformatted')
                 write(iunsterms) ii
                 write(iunsterms) num_nbndv
                 write(iunsterms) n_upper_tot
                 write(iunsterms) numw_prod
                 do jj=1,id_upper
                    write(iunsterms) smat(1:numw_prod, jj)
                 enddo
              else
                 open( unit= iunsterms, file=trim(prefix)//'.vcw_up_overlap.'//nfile, status='old',form='unformatted',&
                      &position='append')
                 do jj=1,id_upper
                    write(iunsterms) smat(1:numw_prod, jj)
                 enddo
              endif
              close(iunsterms)
           endif
           if(id==1) then
              open( unit= iunsterms, file=trim(prefix)//'.vcprim_up.'//nfile, status='unknown',form='unformatted')
              write(iunsterms) ii
              write(iunsterms) num_nbndv
              write(iunsterms) n_upper_tot
              write(iunsterms) numw_prod
              do jj=1,id_upper
                 write(iunsterms) dble(vmat(1:numw_prod, jj))
              enddo
           else
              open( unit= iunsterms, file=trim(prefix)//'.vcprim_up.'//nfile, status='old',form='unformatted',&
                   &position='append')
              do jj=1,id_upper
                 write(iunsterms) dble(vmat(1:numw_prod, jj))
              enddo
           endif
           close(iunsterms)
        endif
        deallocate(smat,vmat)
        deallocate(tmp_sigma)
     enddo
     deallocate(o_mat, o_eigen, o_vector)
  enddo

  deallocate(tmp_g, tmp_prod)
  deallocate(fac)
  deallocate(tmp_reals)
  deallocate(tmp_wp)
  deallocate(tmp_reals_up)
! write global data for number of states and their energies
  iunsterms =  find_free_unit()
  open( unit= iunsterms, file=trim(prefix)//'.upper', status='unknown',form='unformatted')
  write(iunsterms) nbnd
  write(iunsterms) nbnd_normal
  write(iunsterms) num_nbndv
  write(iunsterms) n_upper_tot
  do jj=1,n_upper_tot
     write(iunsterms) ene_reduced(jj)
     write(stdout,*) ene_reduced(jj)
  enddo
  close(iunsterms)
  deallocate(ene_reduced)

! #endif

  return
end subroutine create_upper_states
