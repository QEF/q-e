! FOR GWW
!
! Author: P. Umari
!
 subroutine wannier_valence_terms(nbnd_v,n_set,n_setv)

!this subroutine
!calculates the terms \int dr w^P_i(r)w^P_i(j)*(Psi_v(r)^2)
!put results on wpwp_psi array which is written on disk

   USE io_global,            ONLY : stdout, ionode
   USE io_files,             ONLY : find_free_unit, prefix, iunwfc, nwordwfc,&
                                    diropn, iunigk
   USE mp_global,            ONLY : nproc_pool, me_pool
   USE kinds,                ONLY : DP
   USE basis
   USE klist
   USE constants,            ONLY : e2, pi, tpi, fpi
   USE wvfct,                ONLY : igk, g2kin, npwx, npw, nbnd, nbndx, ecutwfc
   USE control_flags,        ONLY : gamma_only
   USE cell_base,            ONLY : at, alat, tpiba, omega, tpiba2
   USE wannier_gw
   USE fft_base,             ONLY : dffts, dfftp
   USE fft_interfaces,       ONLY : fwfft, invfft
   USE gvect
   USE gvecs,              ONLY : nls, nlsm, doublegrid
   USE uspp
   USE wavefunctions_module, ONLY : psic, evc
   USE realus,               ONLY : adduspos_gamma_r
   USE cell_base,            ONLY : at, bg, omega
   USE mp,                   ONLY : mp_barrier, mp_sum
   USE becmod,               ONLY : calbec

   implicit none


   REAL(kind=DP) :: sca
   INTEGER :: ir




   INTEGER, INTENT(in)  :: nbnd_v !number of KS states considered
   INTEGER, INTENT(in)  :: n_set  !defines the number of states to be read from disk at the same time
   INTEGER, INTENT(in)  :: n_setv  !defines the number of valence states to be read from disk at the same time

   INTEGER :: iungprod, iunuterms

  !  --- Internal definitions ---

   COMPLEX(kind=DP), ALLOCATABLE :: tmpspacei(:,:)
   COMPLEX(kind=DP), ALLOCATABLE :: tmpspacej(:,:)
   COMPLEX(kind=DP), ALLOCATABLE :: tmpspacev(:,:)
   COMPLEX(kind=DP), ALLOCATABLE :: tmpspacec(:)
   INTEGER :: iw,jw, iiw,jjw,jw_begin, hw,hhw
   INTEGER :: ig
   LOGICAL :: exst
   INTEGER :: igk0(npwx)
   REAL(kind=dp) :: g2kin_bp(npwx)
   INTEGER :: npw0
   REAL(kind=DP) :: exxdiv
   INTEGER :: hhv, hv
   REAL(kind=DP), ALLOCATABLE :: wpwp_psi(:,:,:)
   INTEGER :: iunreal, iunterm
   REAL(kind=DP), ALLOCATABLE :: tmpreal1(:),tmpreals1(:),tmpreal2(:),tmpreals2(:)
   REAL(kind=DP), ALLOCATABLE :: becpr(:,:)
   COMPLEX(kind=DP), ALLOCATABLE :: tmpspacei2(:), tmpspacei1(:)


   write(stdout,*) 'Routine wannier_valence_terms : start',nbnd_v,n_set,n_setv


   allocate(wpwp_psi(numw_prod,numw_prod,nbnd_v))
   allocate(tmpreal1(dfftp%nnr),tmpreals1(dffts%nnr))
   allocate(tmpreal2(dfftp%nnr),tmpreals2(dffts%nnr))
   allocate(tmpspacec(dfftp%nnr))
   if(okvan) allocate(becpr(nkb,nbnd_v))

! reads wfcs from iunwfc

   CALL gk_sort(xk(1,1),ngm,g,ecutwfc/tpiba2, &
              &    npw0,igk0,g2kin_bp)
   if(.not.lsmallgrid) then
      allocate(tmpspacei1(ngm),tmpspacei2(ngm))
   else
      allocate(tmpspacei1(npw0),tmpspacei2(npw0))
   endif
   allocate(tmpspacev(npw0,n_setv))


   iungprod = find_free_unit()
   if(.not.lsmallgrid) then
      CALL diropn( iungprod, 'wiwjwfc', ngm*2, exst )
   else
      CALL diropn( iungprod, 'wiwjwfc', npw0*2, exst )
   endif


!sets factors terms

   if(ionode) then
      iunterm =  find_free_unit()
      open( unit= iunterm, file=trim(prefix)//'.wpwp_psi', status='unknown',form='unformatted')
   endif


!open output file

   do hhw=1,ceiling(real(nbnd_v)/real(n_setv))



!      allocate( evc( npwx, nbnd ) )

      write(stdout,*) 'READ HW0',npwx,nbnd,nwordwfc,iunwfc!ATTENZIONE
      call mp_barrier
!      CALL davcio(evc,nwordwfc,iunwfc,1,-1)
      call mp_barrier
       write(stdout,*) 'READ HW1' !ATTENZIONE
       call mp_barrier

      do hw=(hhw-1)*n_setv+1,min(hhw*n_setv,nbnd_v)
          tmpspacev(1:npw0,hw-(hhw-1)*n_setv)=evc(1:npw0,hw)
         if(gstart==2) then
            tmpspacev(1,hw-(hhw-1)*n_setv)=dble(tmpspacev(1,hw-(hhw-1)*n_setv))
         endif
      enddo

      write(stdout,*) 'READ HW'!ATTENZIONE

!      deallocate(evc)

      if  ( nkb > 0 .and. okvan) then
         CALL init_us_2( npw, igk, xk(1,1), vkb )
         call calbec(npw, vkb, tmpspacev, becpr, nbnd_v)
      endif


      do iiw=1,ceiling(real(numw_prod)/real(n_set))
         call mp_barrier
         write(stdout,*) 'READ IW0 ALLOC',iiw!ATTENZIONE
         if(.not.lsmallgrid) then
            allocate(tmpspacei(ngm,n_set))
         else
            allocate(tmpspacei(npw0,n_set))
         endif

         write(stdout,*) 'READ IW0',iiw!ATTENZIONE

         do iw=(iiw-1)*n_set+1,min(iiw*n_set,numw_prod)
            if(.not.lsmallgrid) then
               CALL davcio(tmpspacei(:,iw-(iiw-1)*n_set),ngm*2,iungprod,iw,-1)
            else
               CALL davcio(tmpspacei(:,iw-(iiw-1)*n_set),npw0*2,iungprod,iw,-1)
            endif
            if(gamma_only .and. gstart == 2) then
               tmpspacei(1,iw-(iiw-1)*n_set)=dble(tmpspacei(1,iw-(iiw-1)*n_set))
            endif

            sca=0.d0
            do ir=1,ngm
               sca=sca+2.d0*(tmpspacei(ir,iw-(iiw-1)*n_set)*conjg(tmpspacei(ir,iw-(iiw-1)*n_set)))
            enddo
            if(gstart == 2) sca=sca-tmpspacei(1,iw-(iiw-1)*n_set)*conjg(tmpspacei(1,iw-(iiw-1)*n_set))
            call mp_sum(sca)
            write(stdout,*) 'VERIFICA MOD', iw, iiw, sca


         enddo

         write(stdout,*) 'READ IW'!ATTENZIONE

         do jjw=iiw,ceiling(real(numw_prod)/real(n_set))
            call mp_barrier
            write(stdout,*) 'READ JW0 ALLOC',jjw!ATTENZIONE
            if(.not.lsmallgrid) then
               allocate(tmpspacej(ngm,n_set))
            else
               allocate(tmpspacej(npw0,n_set))
            endif

            write(stdout,*) 'READ JW0',jjw,numw_prod!ATTENZIONE
            do jw=(jjw-1)*n_set+1,min(jjw*n_set,numw_prod)
               if(.not.lsmallgrid) then
                  CALL davcio(tmpspacej(:,jw-(jjw-1)*n_set),ngm*2,iungprod,jw,-1)
               else
                  CALL davcio(tmpspacej(:,jw-(jjw-1)*n_set),npw0*2,iungprod,jw,-1)
               endif
               if(gamma_only .and. gstart == 2) then
                  tmpspacej(1,jw-(jjw-1)*n_set)=dble(tmpspacej(1,jw-(jjw-1)*n_set))
               endif
               call mp_barrier
               write(stdout,*) 'READ JW0',jjw, jw,ngm
            enddo
            call mp_barrier
            write(stdout,*) 'READ JW'!ATTENZIONE

            do hw=(hhw-1)*n_setv+1,min(hhw*n_setv,nbnd_v),2

               psic(:)=(0.d0,0.d0)

               IF ( hw < min(hhw*n_setv,nbnd_v) ) THEN
                  ! ... two ffts at the same time
                  psic(nls(igk(1:npw)))  = tmpspacev(1:npw,hw-(hhw-1)*n_setv) + &
                       ( 0.D0, 1.D0 ) * tmpspacev(1:npw,hw-(hhw-1)*n_setv+1)
                  psic(nlsm(igk(1:npw))) = CONJG( tmpspacev(1:npw,hw-(hhw-1)*n_setv) - &
                       ( 0.D0, 1.D0 ) * tmpspacev(1:npw,hw-(hhw-1)*n_setv+1) )
               ELSE
                  psic(nls(igk(1:npw)))  = tmpspacev(1:npw,hw-(hhw-1)*n_setv)
                  psic(nlsm(igk(1:npw))) = CONJG( tmpspacev(1:npw,hw-(hhw-1)*n_setv) )
               END IF

               CALL invfft ('Wave', psic, dffts)

               tmpreals1(1:dffts%nnr)=dble(psic(1:dffts%nnr))*dble(psic(1:dffts%nnr))
                if(  hw < min(hhw*n_setv,nbnd_v) ) then
                   tmpreals2(1:dffts%nnr)=dimag(psic(1:dffts%nnr))*dimag(psic(1:dffts%nnr))
                endif


                if(doublegrid) then
                   call interpolate(tmpreal1,tmpreals1,1)
                else
                   tmpreal1(:)=tmpreals1(:)
                endif
                if(okvan) call adduspos_gamma_r &
      (hw-(hhw-1)*n_setv,hw-(hhw-1)*n_setv,tmpreal1,1,becpr(:,hw-(hhw-1)*n_setv),becpr(:,hw-(hhw-1)*n_setv))
                if(lsmallgrid) then
                   call interpolate(tmpreal1,tmpreal1,-1)
                endif

                if(  hw < min(hhw*n_setv,nbnd_v) ) then
                   if(doublegrid) then
                      call interpolate(tmpreal2,tmpreals2,1)
                   else
                      tmpreal2(:)=tmpreals2(:)
                   endif
                   if(okvan) call adduspos_gamma_r &
      (hw-(hhw-1)*n_setv+1,hw-(hhw-1)*n_setv+1,tmpreal2,1,becpr(:,hw-(hhw-1)*n_setv+1),becpr(:,hw-(hhw-1)*n_setv+1))
                   if(lsmallgrid) then
                      call interpolate(tmpreal2,tmpreal2,-1)
                   endif
                endif


                do iw=(iiw-1)*n_set+1,min(iiw*n_set,numw_prod)

                   if(.not.lsmallgrid) then

                      tmpspacec(:)=(0.d0,0.d0)
                      tmpspacec(nl(1:ngm))=tmpspacei(1:ngm,iw-(iiw-1)*n_set)
                      tmpspacec(nlm(1:ngm))=CONJG(tmpspacei(1:ngm,iw-(iiw-1)*n_set))

                      CALL invfft ('Dense', tmpspacec, dfftp)

                      tmpspacei1(:)=(0.d0,0.d0)
                      tmpspacec(:)=dble(tmpspacec(:)*tmpreal1(:))
                      CALL fwfft ('Dense', tmpspacec, dfftp)
                      tmpspacei1(1:ngm)=tmpspacec(nl(1:ngm))

                      if(  hw < min(hhw*n_setv,nbnd_v) ) then
                         tmpspacec(:)=(0.d0,0.d0)
                         tmpspacec(nl(1:ngm))=tmpspacei(1:ngm,iw-(iiw-1)*n_set)
                         tmpspacec(nlm(1:ngm))=CONJG(tmpspacei(1:ngm,iw-(iiw-1)*n_set))
                         CALL invfft ('Dense', tmpspacec, dfftp)
                         tmpspacec(:)=dble(tmpspacec(:)*tmpreal2(:))
                         CALL fwfft ('Dense', tmpspacec, dfftp)
                         tmpspacei2(1:ngm)=tmpspacec(nl(1:ngm))
                      endif

                   else
                      tmpspacec(:)=(0.d0,0.d0)
                      tmpspacec(nls(igk0(1:npw0)))=tmpspacei(1:npw0,iw-(iiw-1)*n_set)
                      CALL invfft ('Wave', tmpspacec, dffts)
                      tmpspacec(1:dffts%nnr)=dble(tmpspacec(1:dffts%nnr)*tmpreal1(1:dffts%nnr))
                      CALL fwfft ('Wave', tmpspacec, dffts)
                      tmpspacei1(1:npw0)=tmpspacec(nls(igk0(1:npw0)))
                      if(  hw < min(hhw*n_setv,nbnd_v) ) then
                          tmpspacec(:)=(0.d0,0.d0)
                          tmpspacec(nls(igk0(1:npw0)))=tmpspacei(1:npw0,iw-(iiw-1)*n_set)
                          CALL invfft ('Wave', tmpspacec, dffts)
                          tmpspacec(1:dffts%nnr)=dble(tmpspacec(1:dffts%nnr)*tmpreal2(1:dffts%nnr))
                          CALL fwfft ('Wave', tmpspacec, dffts)
                          tmpspacei2(1:npw0)=tmpspacec(nls(igk0(1:npw0)))
                       endif
                   endif

                   if(iiw==jjw) then
                      jw_begin=iw
                   else
                      jw_begin=(jjw-1)*n_set+1
                   endif
                   do jw=jw_begin,min(jjw*n_set,numw_prod)
                      wpwp_psi(iw,jw,hw)=0.d0
                      if(.not.lsmallgrid) then
                         do ig=1,ngm
                            wpwp_psi(iw,jw,hw)=wpwp_psi(iw,jw,hw)+ &
                                 & 2.d0*dble(conjg(tmpspacei1(ig))*tmpspacej(ig,jw-(jjw-1)*n_set))
                         enddo
                         if(gstart == 2) then
                            wpwp_psi(iw,jw,hw)=wpwp_psi(iw,jw,hw)-&
                                 &dble(conjg(tmpspacei1(1))*tmpspacej(1,jw-(jjw-1)*n_set))
                         endif
                      else
                         do ig=1,npw0
                            wpwp_psi(iw,jw,hw)=wpwp_psi(iw,jw,hw)+ &
                                 & 2.d0*dble(conjg(tmpspacei(ig,iw-(iiw-1)*n_set))*tmpspacej(ig,jw-(jjw-1)*n_set))
                         enddo
                         if(gstart == 2) then
                            wpwp_psi(iw,jw,hw)=wpwp_psi(iw,jw,hw)-&
                                 &dble(conjg(tmpspacei(1,iw-(iiw-1)*n_set))*tmpspacej(1,jw-(jjw-1)*n_set))
                         endif
                      endif
                      call mp_sum(wpwp_psi(iw,jw,hw))

                      wpwp_psi(jw,iw,hw)=wpwp_psi(iw,jw,hw)
                      if(  hw < min(hhw*n_setv,nbnd_v) ) then
                         wpwp_psi(iw,jw,hw+1)=0.d0
                         if(.not.lsmallgrid) then
                            do ig=1,ngm
                               wpwp_psi(iw,jw,hw+1)=wpwp_psi(iw,jw,hw+1)+ &
                                    & 2.d0*dble(conjg(tmpspacei2(ig))*tmpspacej(ig,jw-(jjw-1)*n_set))
                            enddo
                            if(gstart == 2) then
                               wpwp_psi(iw,jw,hw+1)=wpwp_psi(iw,jw,hw+1)-&
                                    &dble(conjg(tmpspacei2(1))*tmpspacej(1,jw-(jjw-1)*n_set))
                            endif
                         else
                            do ig=1,npw0
                               wpwp_psi(iw,jw,hw+1)=wpwp_psi(iw,jw,hw+1)+ &
                                    & 2.d0*dble(conjg(tmpspacei2(ig))*tmpspacej(ig,jw-(jjw-1)*n_set))
                            enddo
                            if(gstart == 2) then
                               wpwp_psi(iw,jw,hw+1)=wpwp_psi(iw,jw,hw+1)-&
                                    &dble(conjg(tmpspacei2(1))*tmpspacej(1,jw-(jjw-1)*n_set))
                            endif
                         endif
                         call mp_sum(wpwp_psi(iw,jw,hw+1))
                         wpwp_psi(jw,iw,hw+1)=wpwp_psi(iw,jw,hw+1)
                      endif
                  enddo
                enddo
            enddo
            deallocate(tmpspacej)
         enddo
         deallocate(tmpspacei)
      enddo
   enddo
   if(ionode) then
      write(iunterm) numw_prod
      write(iunterm) nbnd_v
      do hw=1,nbnd_v
         do iw=1,numw_prod
            write(iunterm) wpwp_psi(iw,1:iw,hw)
         enddo
      enddo
   endif
   close(iungprod)
   if(ionode) close(iunterm)
   deallocate(tmpspacev,wpwp_psi)
   deallocate(tmpreal1, tmpreals1)
   deallocate(tmpreal2, tmpreals2)
   deallocate(tmpspacec)
   deallocate(tmpspacei1,tmpspacei2)
   return
 end subroutine wannier_valence_terms



 subroutine wannier_valence_terms_cutoff(nbnd_v,n_set,n_setv)

!this subroutine
!calculates the terms \int dr w^P_i(r)w^P_i(j)*(Psi_v(r)^2)
!put results on wpwp_psi array which is written on disk

   USE io_global,            ONLY : stdout, ionode
   USE io_files,             ONLY : find_free_unit, prefix, iunwfc, nwordwfc, iunigk
   USE mp_global,            ONLY : nproc_pool, me_pool
   USE kinds,                ONLY : DP
   USE basis
   USE klist
   USE constants,            ONLY : e2, pi, tpi, fpi
   USE wvfct,                ONLY : igk, g2kin, npwx, npw, nbnd, nbndx, ecutwfc
   USE control_flags,        ONLY : gamma_only
   USE cell_base,            ONLY : at, alat, tpiba, omega, tpiba2
   USE wannier_gw
   USE fft_base,             ONLY : dffts, dfftp
   USE fft_interfaces,       ONLY : fwfft, invfft
   USE gvect
   USE gvecs,              ONLY : nls, nlsm, doublegrid
   USE uspp
   USE wavefunctions_module, ONLY : psic, evc
   USE realus,               ONLY : adduspos_gamma_r
   USE cell_base,            ONLY : at, bg, omega
   USE mp,                   ONLY : mp_barrier, mp_sum
   USE becmod,               ONLY : calbec

   implicit none




   REAL(kind=DP) :: sca
   INTEGER :: ir




   INTEGER, INTENT(in)  :: nbnd_v !number of KS states considered
   INTEGER, INTENT(in)  :: n_set  !defines the number of states to be read from disk at the same time
   INTEGER, INTENT(in)  :: n_setv  !defines the number of valence states to be read from disk at the same time

   INTEGER :: iungprod, iunuterms

  !  --- Internal definitions ---

   COMPLEX(kind=DP), ALLOCATABLE :: tmpspacei(:,:)
   COMPLEX(kind=DP), ALLOCATABLE :: tmpspacej(:,:)
   COMPLEX(kind=DP), ALLOCATABLE :: tmpspacev(:,:)
   COMPLEX(kind=DP), ALLOCATABLE :: tmpspaced(:,:)
   COMPLEX(kind=DP), ALLOCATABLE :: tmpspacec(:)
   INTEGER :: iw,jw, iiw,jjw,jw_begin, hw,hhw, iw_begin
   INTEGER :: ig
   LOGICAL :: exst
   INTEGER :: igk0(npwx)
   REAL(kind=dp) :: g2kin_bp(npwx)
   INTEGER :: npw0
   REAL(kind=DP) :: exxdiv
   INTEGER :: hhv, hv
   REAL(kind=DP) :: wpwp_psi
   INTEGER :: iunreal, iunterm, iunterm2
   REAL(kind=DP), ALLOCATABLE :: tmpreal1(:),tmpreals1(:),tmpreal2(:),tmpreals2(:)
   REAL(kind=DP), ALLOCATABLE :: becpr(:,:)
   COMPLEX(kind=DP), ALLOCATABLE :: tmpspacei2(:), tmpspacei1(:)
   INTEGER, ALLOCATABLE :: ok_table(:,:)
   INTEGER, ALLOCATABLE :: index_terms(:,:)

   write(stdout,*) 'Routine wannier_valence_terms : start',nbnd_v,n_set,n_setv

   allocate(ok_table(n_set, n_set))

   allocate(tmpreal1(dfftp%nnr),tmpreals1(dffts%nnr))
   allocate(tmpreal2(dfftp%nnr),tmpreals2(dffts%nnr))
   allocate(tmpspaced(dfftp%nnr,n_set))
   allocate(tmpspacec(dfftp%nnr))


   numw_prodprod=0

   if(okvan) allocate(becpr(nkb,nbnd_v))

! reads wfcs from iunwfc

   CALL gk_sort(xk(1,1),ngm,g,ecutwfc/tpiba2, &
              &    npw0,igk0,g2kin_bp)
   if(.not.lsmallgrid) then
      allocate(tmpspacei1(ngm),tmpspacei2(ngm))
   else
      allocate(tmpspacei1(npw0),tmpspacei2(npw0))
   endif
   allocate(tmpspacev(npw0,n_setv))


   iungprod = find_free_unit()
   if(.not.lsmallgrid) then
      CALL diropn( iungprod, 'wiwjwfc', ngm*2, exst )
   else
      CALL diropn( iungprod, 'wiwjwfc', npw0*2, exst )
   endif


!sets factors terms

   if(ionode) then
      iunterm =  find_free_unit()
      open( unit= iunterm, file=trim(prefix)//'.wpwp_psi', status='unknown',form='unformatted')
      iunterm2 = find_free_unit()
      open( unit= iunterm2, file=trim(prefix)//'.wpwp_psi_index', status='unknown',form='unformatted')
      write(iunterm2) numw_prod
      write(iunterm2) nbnd_v
      write(iunterm2) numw_prod!at the end this must be re-written as numw_prodprod
   endif



!open output file

   do hhw=1,ceiling(real(nbnd_v)/real(n_setv))



!      allocate( evc( npwx, nbnd ) )

      write(stdout,*) 'READ HW0',npwx,nbnd,nwordwfc,iunwfc!ATTENZIONE
      call mp_barrier
!      CALL davcio(evc,nwordwfc,iunwfc,1,-1)
      call mp_barrier
       write(stdout,*) 'READ HW1' !ATTENZIONE
       call mp_barrier

      do hw=(hhw-1)*n_setv+1,min(hhw*n_setv,nbnd_v)
          tmpspacev(1:npw0,hw-(hhw-1)*n_setv)=evc(1:npw0,hw)
         if(gstart==2) then
            tmpspacev(1,hw-(hhw-1)*n_setv)=dble(tmpspacev(1,hw-(hhw-1)*n_setv))
         endif
      enddo

      write(stdout,*) 'READ HW'!ATTENZIONE

!      deallocate(evc)

      if  ( nkb > 0 .and. okvan) then
         CALL init_us_2( npw, igk, xk(1,1), vkb )
         call calbec(npw, vkb, tmpspacev, becpr, nbnd_v)
      endif


      do iiw=1,ceiling(real(numw_prod)/real(n_set))
         call mp_barrier
         write(stdout,*) 'READ IW0 ALLOC',iiw!ATTENZIONE
         if(.not.lsmallgrid) then
            allocate(tmpspacei(ngm,n_set))
         else
            allocate(tmpspacei(npw0,n_set))
         endif

         write(stdout,*) 'READ IW0',iiw!ATTENZIONE

         do iw=(iiw-1)*n_set+1,min(iiw*n_set,numw_prod)
            if(.not.lsmallgrid) then
               CALL davcio(tmpspacei(:,iw-(iiw-1)*n_set),ngm*2,iungprod,iw,-1)
            else
               CALL davcio(tmpspacei(:,iw-(iiw-1)*n_set),npw0*2,iungprod,iw,-1)
            endif
            if(gamma_only .and. gstart == 2) then
               tmpspacei(1,iw-(iiw-1)*n_set)=dble(tmpspacei(1,iw-(iiw-1)*n_set))
            endif

            sca=0.d0
            do ir=1,ngm
               sca=sca+2.d0*(tmpspacei(ir,iw-(iiw-1)*n_set)*conjg(tmpspacei(ir,iw-(iiw-1)*n_set)))
            enddo
            if(gstart == 2) sca=sca-tmpspacei(1,iw-(iiw-1)*n_set)*conjg(tmpspacei(1,iw-(iiw-1)*n_set))
            call mp_sum(sca)
            write(stdout,*) 'VERIFICA MOD', iw, iiw, sca


         enddo

         write(stdout,*) 'READ IW'!ATTENZIONE

         do jjw=iiw,ceiling(real(numw_prod)/real(n_set))
            call mp_barrier
            write(stdout,*) 'READ JW0 ALLOC',jjw!ATTENZIONE
            if(.not.lsmallgrid) then
               allocate(tmpspacej(ngm,n_set))
            else
               allocate(tmpspacej(npw0,n_set))
            endif

            write(stdout,*) 'READ JW0',jjw,numw_prod!ATTENZIONE
            do jw=(jjw-1)*n_set+1,min(jjw*n_set,numw_prod)
               if(.not.lsmallgrid) then
                  CALL davcio(tmpspacej(:,jw-(jjw-1)*n_set),ngm*2,iungprod,jw,-1)
               else
                  CALL davcio(tmpspacej(:,jw-(jjw-1)*n_set),npw0*2,iungprod,jw,-1)
               endif
               if(gamma_only .and. gstart == 2) then
                  tmpspacej(1,jw-(jjw-1)*n_set)=dble(tmpspacej(1,jw-(jjw-1)*n_set))
               endif
               call mp_barrier
               write(stdout,*) 'READ JW0',jjw, jw,ngm
            enddo
            call mp_barrier
            write(stdout,*) 'READ JW'!ATTENZIONE

!look for non negligeable components
            ok_table(:,:)=0
            do iw=(iiw-1)*n_set+1,min(iiw*n_set,numw_prod)

               if(.not.lsmallgrid) then

                  tmpspaced(:,iw-(iiw-1)*n_set)=(0.d0,0.d0)
                  tmpspaced(nl(1:ngm), iw-(iiw-1)*n_set)=tmpspacei(1:ngm,iw-(iiw-1)*n_set)
                  tmpspaced(nlm(1:ngm),iw-(iiw-1)*n_set)=CONJG(tmpspacei(1:ngm,iw-(iiw-1)*n_set))

                  CALL invfft ('Dense', tmpspaced(:,iw-(iiw-1)*n_set), dfftp)
               else
                  write(stdout,*) 'lsmallgrid not implemented'
                  stop
               endif
            end do

            do jw=(jjw-1)*n_set+1,min(jjw*n_set,numw_prod)
               if(.not.lsmallgrid) then
                  tmpspacec(:)=(0.d0,0.d0)
                  tmpspacec(nl(1:ngm))=tmpspacej(1:ngm,jw-(jjw-1)*n_set)
                  tmpspacec(nlm(1:ngm))=CONJG(tmpspacej(1:ngm,jw-(jjw-1)*n_set))
                  CALL invfft ('Dense', tmpspacec, dfftp)
               else
                  write(stdout,*) 'lsmallgrid not implemented'
                  stop
               endif
               if(iiw==jjw) then
                  iw_begin=jw
               else
                  iw_begin=(iiw-1)*n_set+1
               endif
               do iw=iw_begin,min(iiw*n_set,numw_prod)

                  sca=0.d0
                  do ir=1,dfftp%nnr
                     sca=sca+(tmpspacec(ir)*tmpspaced(ir, iw-(iiw-1)*n_set))**2.d0
                  enddo
                  call mp_sum(sca)
                  if(sca >= cutoff_wpr_wpr) then
                     !if(ionode) write(stdout,*) 'PASSED WPWP', numw_prodprod+1,sca
                     numw_prodprod=numw_prodprod+1
                     ok_table(iw-(iiw-1)*n_set, jw-(jjw-1)*n_set)=numw_prodprod
                     if(iiw==jjw) then
                        ok_table(jw-(jjw-1)*n_set, iw-(iiw-1)*n_set)=numw_prodprod
                     endif
                     if(ionode) write(iunterm2) iw,jw
                  endif
               enddo
            enddo




            do hw=(hhw-1)*n_setv+1,min(hhw*n_setv,nbnd_v),2

               psic(:)=(0.d0,0.d0)

               IF ( hw < min(hhw*n_setv,nbnd_v) ) THEN
                  ! ... two ffts at the same time
                  psic(nls(igk(1:npw)))  = tmpspacev(1:npw,hw-(hhw-1)*n_setv) + &
                       ( 0.D0, 1.D0 ) * tmpspacev(1:npw,hw-(hhw-1)*n_setv+1)
                  psic(nlsm(igk(1:npw))) = CONJG( tmpspacev(1:npw,hw-(hhw-1)*n_setv) - &
                       ( 0.D0, 1.D0 ) * tmpspacev(1:npw,hw-(hhw-1)*n_setv+1) )
               ELSE
                  psic(nls(igk(1:npw)))  = tmpspacev(1:npw,hw-(hhw-1)*n_setv)
                  psic(nlsm(igk(1:npw))) = CONJG( tmpspacev(1:npw,hw-(hhw-1)*n_setv) )
               END IF

               CALL invfft ('Wave', psic, dffts)

               tmpreals1(1:dffts%nnr)=dble(psic(1:dffts%nnr))*dble(psic(1:dffts%nnr))
                if(  hw < min(hhw*n_setv,nbnd_v) ) then
                   tmpreals2(1:dffts%nnr)=dimag(psic(1:dffts%nnr))*dimag(psic(1:dffts%nnr))
                endif


                if(doublegrid) then
                   call interpolate(tmpreal1,tmpreals1,1)
                else
                   tmpreal1(:)=tmpreals1(:)
                endif
                if(okvan) call adduspos_gamma_r &
      (hw-(hhw-1)*n_setv,hw-(hhw-1)*n_setv,tmpreal1,1,becpr(:,hw-(hhw-1)*n_setv),becpr(:,hw-(hhw-1)*n_setv))
                if(lsmallgrid) then
                   call interpolate(tmpreal1,tmpreal1,-1)
                endif

                if(  hw < min(hhw*n_setv,nbnd_v) ) then
                   if(doublegrid) then
                      call interpolate(tmpreal2,tmpreals2,1)
                   else
                      tmpreal2(:)=tmpreals2(:)
                   endif
                   if(okvan) call adduspos_gamma_r &
      (hw-(hhw-1)*n_setv+1,hw-(hhw-1)*n_setv+1,tmpreal2,1,becpr(:,hw-(hhw-1)*n_setv+1),becpr(:,hw-(hhw-1)*n_setv+1))
                   if(lsmallgrid) then
                      call interpolate(tmpreal2,tmpreal2,-1)
                   endif
                endif


                do iw=(iiw-1)*n_set+1,min(iiw*n_set,numw_prod)

                   if(.not.lsmallgrid) then

                      tmpspacec(:)=(0.d0,0.d0)
                      tmpspacec(nl(1:ngm))=tmpspacei(1:ngm,iw-(iiw-1)*n_set)
                      tmpspacec(nlm(1:ngm))=CONJG(tmpspacei(1:ngm,iw-(iiw-1)*n_set))

                      CALL invfft ('Dense', tmpspacec, dfftp)

                      tmpspacei1(:)=(0.d0,0.d0)
                      tmpspacec(:)=dble(tmpspacec(:)*tmpreal1(:))
                      CALL fwfft ('Dense', tmpspacec, dfftp)
                      tmpspacei1(1:ngm)=tmpspacec(nl(1:ngm))

                      if(  hw < min(hhw*n_setv,nbnd_v) ) then
                         tmpspacec(:)=(0.d0,0.d0)
                         tmpspacec(nl(1:ngm))=tmpspacei(1:ngm,iw-(iiw-1)*n_set)
                         tmpspacec(nlm(1:ngm))=CONJG(tmpspacei(1:ngm,iw-(iiw-1)*n_set))
                         CALL invfft ('Dense', tmpspacec, dfftp)
                         tmpspacec(:)=dble(tmpspacec(:)*tmpreal2(:))
                         CALL fwfft ('Dense', tmpspacec, dfftp)
                         tmpspacei2(1:ngm)=tmpspacec(nl(1:ngm))
                      endif

                   else
                      tmpspacec(:)=(0.d0,0.d0)
                      tmpspacec(nls(igk0(1:npw0)))=tmpspacei(1:npw0,iw-(iiw-1)*n_set)
                      CALL invfft ('Wave', tmpspacec, dffts)
                      tmpspacec(1:dffts%nnr)=dble(tmpspacec(1:dffts%nnr)*tmpreal1(1:dffts%nnr))
                      CALL fwfft ('Wave', tmpspacec, dffts)
                      tmpspacei1(1:npw0)=tmpspacec(nls(igk0(1:npw0)))
                      if(  hw < min(hhw*n_setv,nbnd_v) ) then
                          tmpspacec(:)=(0.d0,0.d0)
                          tmpspacec(nls(igk0(1:npw0)))=tmpspacei(1:npw0,iw-(iiw-1)*n_set)
                          CALL invfft ('Wave', tmpspacec, dffts)
                          tmpspacec(1:dffts%nnr)=dble(tmpspacec(1:dffts%nnr)*tmpreal2(1:dffts%nnr))
                          CALL fwfft ('Wave', tmpspacec, dffts)
                          tmpspacei2(1:npw0)=tmpspacec(nls(igk0(1:npw0)))
                       endif
                   endif

                   if(iiw==jjw) then
                      jw_begin=iw
                   else
                      jw_begin=(jjw-1)*n_set+1
                   endif
                   do jw=jw_begin,min(jjw*n_set,numw_prod)
                      if(ok_table(iw-(iiw-1)*n_set,jw-(jjw-1)*n_set) /= 0 ) then
                         wpwp_psi=0.d0
                         if(.not.lsmallgrid) then
                            do ig=1,ngm
                               wpwp_psi=wpwp_psi+ &
                                    & 2.d0*dble(conjg(tmpspacei1(ig))*tmpspacej(ig,jw-(jjw-1)*n_set))
                            enddo
                            if(gstart == 2) then
                               wpwp_psi=wpwp_psi-&
                                    &dble(conjg(tmpspacei1(1))*tmpspacej(1,jw-(jjw-1)*n_set))
                            endif
                         else
                            do ig=1,npw0
                               wpwp_psi=wpwp_psi+ &
                                    & 2.d0*dble(conjg(tmpspacei(ig,iw-(iiw-1)*n_set))*tmpspacej(ig,jw-(jjw-1)*n_set))
                            enddo
                            if(gstart == 2) then
                               wpwp_psi=wpwp_psi-&
                                    &dble(conjg(tmpspacei(1,iw-(iiw-1)*n_set))*tmpspacej(1,jw-(jjw-1)*n_set))
                            endif
                         endif
                         call mp_sum(wpwp_psi)
                         if(ionode) write(iunterm) ok_table(iw-(iiw-1)*n_set,jw-(jjw-1)*n_set),hw,wpwp_psi
                         !SCRIVI wpwp_psi su file
                         if(  hw < min(hhw*n_setv,nbnd_v) ) then
                            wpwp_psi=0.d0
                            if(.not.lsmallgrid) then
                               do ig=1,ngm
                                  wpwp_psi=wpwp_psi+ &
                                       & 2.d0*dble(conjg(tmpspacei2(ig))*tmpspacej(ig,jw-(jjw-1)*n_set))
                               enddo
                               if(gstart == 2) then
                                  wpwp_psi=wpwp_psi-&
                                       &dble(conjg(tmpspacei2(1))*tmpspacej(1,jw-(jjw-1)*n_set))
                               endif
                            else
                               do ig=1,npw0
                                  wpwp_psi=wpwp_psi+ &
                                       & 2.d0*dble(conjg(tmpspacei2(ig))*tmpspacej(ig,jw-(jjw-1)*n_set))
                               enddo
                               if(gstart == 2) then
                                  wpwp_psi=wpwp_psi-&
                                       &dble(conjg(tmpspacei2(1))*tmpspacej(1,jw-(jjw-1)*n_set))
                               endif
                            endif
                            call mp_sum(wpwp_psi)
                            if(ionode) write(iunterm) ok_table(iw-(iiw-1)*n_set,jw-(jjw-1)*n_set),hw+1,wpwp_psi
                         endif
                      endif
                   enddo
                enddo
             enddo
             deallocate(tmpspacej)
          enddo
          deallocate(tmpspacei)
       enddo
    enddo
    close(iungprod)
    if(ionode) then
       close(iunterm)
       close(iunterm2)
       allocate(index_terms(2,numw_prodprod))
       open( unit= iunterm2, file=trim(prefix)//'.wpwp_psi_index', status='old',action='readwrite',form='unformatted')
       rewind(iunterm2)
       read(iunterm2) iw!
       read(iunterm2) iw!nbnd_v
       read(iunterm2) iw!numw_prodprod
       do ig=1, numw_prodprod
          read(iunterm2) index_terms(1,ig), index_terms(2,ig)
       enddo
       rewind(iunterm2)
       write(iunterm2) numw_prod
       write(iunterm2) nbnd_v
       write(iunterm2) numw_prodprod
       do ig=1, numw_prodprod
          write(iunterm2) index_terms(1,ig), index_terms(2,ig)
       enddo
       write(stdout,*) 'TOTAL NUMBER OF PRODUCTS OF WANNIER PRODUCTS :', numw_prodprod
       close(iunterm2)
       deallocate(index_terms)
    endif
    deallocate(tmpspacev)
    deallocate(tmpreal1, tmpreals1)
    deallocate(tmpreal2, tmpreals2)
    deallocate(tmpspacec)
    deallocate(tmpspacei1,tmpspacei2)
    deallocate(tmpspaced)
    return
  end subroutine wannier_valence_terms_cutoff
 !
 !
 subroutine wannier_valence_terms_distance(nbnd_v,n_set,n_setv)

!this subroutine
!calculates the terms \int dr w^P_i(r)w^P_i(j)*(Psi_v(r)^2)
!put results on wpwp_psi array which is written on disk
!cutoff at distance

   USE io_global,            ONLY : stdout, ionode
   USE io_files,             ONLY : find_free_unit, prefix, iunwfc, nwordwfc,&
                                    diropn, iunigk
   USE mp_global,            ONLY : nproc_pool, me_pool
   USE kinds,                ONLY : DP
   USE basis
   USE klist
   USE constants,            ONLY : e2, pi, tpi, fpi
   USE wvfct,                ONLY : igk, g2kin, npwx, npw, nbnd, nbndx, ecutwfc
   USE control_flags,        ONLY : gamma_only
   USE cell_base,            ONLY : at, alat, tpiba, omega, tpiba2
   USE wannier_gw
   USE fft_base,             ONLY : dffts, dfftp
   USE fft_interfaces,       ONLY : fwfft, invfft
   USE gvect
   USE gvecs,              ONLY : nls, nlsm, doublegrid
   USE uspp
   USE wavefunctions_module, ONLY : psic, evc
   USE realus,               ONLY : adduspos_gamma_r
   USE cell_base,            ONLY : at, bg, omega
   USE mp,                   ONLY : mp_barrier, mp_sum
   USE becmod,               ONLY : calbec

   implicit none

   REAL(kind=DP) :: sca
   INTEGER :: ir




   INTEGER, INTENT(in)  :: nbnd_v !number of KS states considered
   INTEGER, INTENT(in)  :: n_set  !defines the number of states to be read from disk at the same time
   INTEGER, INTENT(in)  :: n_setv  !defines the number of valence states to be read from disk at the same time

   INTEGER :: iungprod, iunuterms

  !  --- Internal definitions ---

   COMPLEX(kind=DP), ALLOCATABLE :: tmpspacei(:,:)
   COMPLEX(kind=DP), ALLOCATABLE :: tmpspacej(:,:)
   COMPLEX(kind=DP), ALLOCATABLE :: tmpspacev(:,:)
   COMPLEX(kind=DP), ALLOCATABLE :: tmpspaced(:,:)
   COMPLEX(kind=DP), ALLOCATABLE :: tmpspacec(:)
   INTEGER :: iw,jw, iiw,jjw,jw_begin, hw,hhw, iw_begin
   INTEGER :: ig
   LOGICAL :: exst
   INTEGER :: igk0(npwx)
   REAL(kind=dp) :: g2kin_bp(npwx)
   INTEGER :: npw0
   REAL(kind=DP) :: exxdiv
   INTEGER :: hhv, hv
   REAL(kind=DP) :: wpwp_psi
   INTEGER :: iunreal, iunterm, iunterm2
   REAL(kind=DP), ALLOCATABLE :: tmpreal1(:),tmpreals1(:),tmpreal2(:),tmpreals2(:)
   REAL(kind=DP), ALLOCATABLE :: becpr(:,:)
   COMPLEX(kind=DP), ALLOCATABLE :: tmpspacei2(:), tmpspacei1(:)
   INTEGER, ALLOCATABLE :: ok_table(:,:)
   INTEGER, ALLOCATABLE :: index_terms(:,:)

   write(stdout,*) 'Routine wannier_valence_terms : start',nbnd_v,n_set,n_setv

   allocate(ok_table(n_set, n_set))

   allocate(tmpreal1(dfftp%nnr),tmpreals1(dffts%nnr))
   allocate(tmpreal2(dfftp%nnr),tmpreals2(dffts%nnr))
   allocate(tmpspaced(dfftp%nnr,n_set))
   allocate(tmpspacec(dfftp%nnr))


   numw_prodprod=0

   if(okvan) allocate(becpr(nkb,nbnd_v))

! reads wfcs from iunwfc

   CALL gk_sort(xk(1,1),ngm,g,ecutwfc/tpiba2, &
              &    npw0,igk0,g2kin_bp)
   if(.not.lsmallgrid) then
      allocate(tmpspacei1(ngm),tmpspacei2(ngm))
   else
      allocate(tmpspacei1(npw0),tmpspacei2(npw0))
   endif
   allocate(tmpspacev(npw0,n_setv))


   iungprod = find_free_unit()
   if(.not.lsmallgrid) then
      CALL diropn( iungprod, 'wiwjwfc', ngm*2, exst )
   else
      CALL diropn( iungprod, 'wiwjwfc', npw0*2, exst )
   endif


!sets factors terms

   if(ionode) then
      iunterm =  find_free_unit()
      open( unit= iunterm, file=trim(prefix)//'.wpwp_psi', status='unknown',form='unformatted')
      iunterm2 = find_free_unit()
      open( unit= iunterm2, file=trim(prefix)//'.wpwp_psi_index', status='unknown',form='unformatted')
      write(iunterm2) numw_prod
      write(iunterm2) nbnd_v
      write(iunterm2) numw_prod!at the end this must be re-written as numw_prodprod
   endif



!open output file

   do hhw=1,ceiling(real(nbnd_v)/real(n_setv))



!      allocate( evc( npwx, nbnd ) )

      write(stdout,*) 'READ HW0',npwx,nbnd,nwordwfc,iunwfc!ATTENZIONE
      call mp_barrier
!      CALL davcio(evc,nwordwfc,iunwfc,1,-1)
      call mp_barrier
       write(stdout,*) 'READ HW1' !ATTENZIONE
       call mp_barrier

      do hw=(hhw-1)*n_setv+1,min(hhw*n_setv,nbnd_v)
          tmpspacev(1:npw0,hw-(hhw-1)*n_setv)=evc(1:npw0,hw)
         if(gstart==2) then
            tmpspacev(1,hw-(hhw-1)*n_setv)=dble(tmpspacev(1,hw-(hhw-1)*n_setv))
         endif
      enddo

      write(stdout,*) 'READ HW'!ATTENZIONE

!      deallocate(evc)

      if  ( nkb > 0 .and. okvan) then
         CALL init_us_2( npw, igk, xk(1,1), vkb )
         call calbec(npw, vkb, tmpspacev, becpr, nbnd_v)
      endif


      do iiw=1,ceiling(real(numw_prod)/real(n_set))
         call mp_barrier
         write(stdout,*) 'READ IW0 ALLOC',iiw!ATTENZIONE
         if(.not.lsmallgrid) then
            allocate(tmpspacei(ngm,n_set))
         else
            allocate(tmpspacei(npw0,n_set))
         endif

         write(stdout,*) 'READ IW0',iiw!ATTENZIONE

         do iw=(iiw-1)*n_set+1,min(iiw*n_set,numw_prod)
            if(.not.lsmallgrid) then
               CALL davcio(tmpspacei(:,iw-(iiw-1)*n_set),ngm*2,iungprod,iw,-1)
            else
               CALL davcio(tmpspacei(:,iw-(iiw-1)*n_set),npw0*2,iungprod,iw,-1)
            endif
            if(gamma_only .and. gstart == 2) then
               tmpspacei(1,iw-(iiw-1)*n_set)=dble(tmpspacei(1,iw-(iiw-1)*n_set))
            endif

            sca=0.d0
            do ir=1,ngm
               sca=sca+2.d0*(tmpspacei(ir,iw-(iiw-1)*n_set)*conjg(tmpspacei(ir,iw-(iiw-1)*n_set)))
            enddo
            if(gstart == 2) sca=sca-tmpspacei(1,iw-(iiw-1)*n_set)*conjg(tmpspacei(1,iw-(iiw-1)*n_set))
            call mp_sum(sca)
            write(stdout,*) 'VERIFICA MOD', iw, iiw, sca


         enddo

         write(stdout,*) 'READ IW'!ATTENZIONE

         do jjw=iiw,ceiling(real(numw_prod)/real(n_set))
            call mp_barrier
            write(stdout,*) 'READ JW0 ALLOC',jjw!ATTENZIONE
            if(.not.lsmallgrid) then
               allocate(tmpspacej(ngm,n_set))
            else
               allocate(tmpspacej(npw0,n_set))
            endif

            write(stdout,*) 'READ JW0',jjw,numw_prod!ATTENZIONE
            do jw=(jjw-1)*n_set+1,min(jjw*n_set,numw_prod)
               if(.not.lsmallgrid) then
                  CALL davcio(tmpspacej(:,jw-(jjw-1)*n_set),ngm*2,iungprod,jw,-1)
               else
                  CALL davcio(tmpspacej(:,jw-(jjw-1)*n_set),npw0*2,iungprod,jw,-1)
               endif
               if(gamma_only .and. gstart == 2) then
                  tmpspacej(1,jw-(jjw-1)*n_set)=dble(tmpspacej(1,jw-(jjw-1)*n_set))
               endif
               call mp_barrier
               write(stdout,*) 'READ JW0',jjw, jw,ngm
            enddo
            call mp_barrier
            write(stdout,*) 'READ JW'!ATTENZIONE

!look for non negligeable components
            ok_table(:,:)=0
            do iw=(iiw-1)*n_set+1,min(iiw*n_set,numw_prod)

               if(.not.lsmallgrid) then

                  tmpspaced(:,iw-(iiw-1)*n_set)=(0.d0,0.d0)
                  tmpspaced(nl(1:ngm), iw-(iiw-1)*n_set)=tmpspacei(1:ngm,iw-(iiw-1)*n_set)
                  tmpspaced(nlm(1:ngm),iw-(iiw-1)*n_set)=CONJG(tmpspacei(1:ngm,iw-(iiw-1)*n_set))

                  CALL invfft ('Dense', tmpspaced(:,iw-(iiw-1)*n_set), dfftp)
               else
                  write(stdout,*) 'lsmallgrid not implemented'
                  stop
               endif
            end do

            do jw=(jjw-1)*n_set+1,min(jjw*n_set,numw_prod)
               if(.not.lsmallgrid) then
                  tmpspacec(:)=(0.d0,0.d0)
                  tmpspacec(nl(1:ngm))=tmpspacej(1:ngm,jw-(jjw-1)*n_set)
                  tmpspacec(nlm(1:ngm))=CONJG(tmpspacej(1:ngm,jw-(jjw-1)*n_set))
                  CALL invfft ('Dense', tmpspacec, dfftp)
               else
                  write(stdout,*) 'lsmallgrid not implemented'
                  stop
               endif
               if(iiw==jjw) then
                  iw_begin=jw
               else
                  iw_begin=(iiw-1)*n_set+1
               endif
               do iw=iw_begin,min(iiw*n_set,numw_prod)

                  if(l_on_products(iw,jw)) then
                     numw_prodprod=numw_prodprod+1
                     ok_table(iw-(iiw-1)*n_set, jw-(jjw-1)*n_set)=numw_prodprod
                     if(iiw==jjw) then
                        ok_table(jw-(jjw-1)*n_set, iw-(iiw-1)*n_set)=numw_prodprod
                     endif
                     if(ionode) write(iunterm2) iw,jw
                  endif
               enddo
            enddo




            do hw=(hhw-1)*n_setv+1,min(hhw*n_setv,nbnd_v),2

               psic(:)=(0.d0,0.d0)

               IF ( hw < min(hhw*n_setv,nbnd_v) ) THEN
                  ! ... two ffts at the same time
                  psic(nls(igk(1:npw)))  = tmpspacev(1:npw,hw-(hhw-1)*n_setv) + &
                       ( 0.D0, 1.D0 ) * tmpspacev(1:npw,hw-(hhw-1)*n_setv+1)
                  psic(nlsm(igk(1:npw))) = CONJG( tmpspacev(1:npw,hw-(hhw-1)*n_setv) - &
                       ( 0.D0, 1.D0 ) * tmpspacev(1:npw,hw-(hhw-1)*n_setv+1) )
               ELSE
                  psic(nls(igk(1:npw)))  = tmpspacev(1:npw,hw-(hhw-1)*n_setv)
                  psic(nlsm(igk(1:npw))) = CONJG( tmpspacev(1:npw,hw-(hhw-1)*n_setv) )
               END IF

               CALL invfft ('Wave', psic, dffts)

               tmpreals1(1:dffts%nnr)=dble(psic(1:dffts%nnr))*dble(psic(1:dffts%nnr))
                if(  hw < min(hhw*n_setv,nbnd_v) ) then
                   tmpreals2(1:dffts%nnr)=dimag(psic(1:dffts%nnr))*dimag(psic(1:dffts%nnr))
                endif


                if(doublegrid) then
                   call interpolate(tmpreal1,tmpreals1,1)
                else
                   tmpreal1(:)=tmpreals1(:)
                endif
                if(okvan) call adduspos_gamma_r &
      (hw-(hhw-1)*n_setv,hw-(hhw-1)*n_setv,tmpreal1,1,becpr(:,hw-(hhw-1)*n_setv),becpr(:,hw-(hhw-1)*n_setv))
                if(lsmallgrid) then
                   call interpolate(tmpreal1,tmpreal1,-1)
                endif

                if(  hw < min(hhw*n_setv,nbnd_v) ) then
                   if(doublegrid) then
                      call interpolate(tmpreal2,tmpreals2,1)
                   else
                      tmpreal2(:)=tmpreals2(:)
                   endif
                   if(okvan) call adduspos_gamma_r &
      (hw-(hhw-1)*n_setv+1,hw-(hhw-1)*n_setv+1,tmpreal2,1,becpr(:,hw-(hhw-1)*n_setv+1),becpr(:,hw-(hhw-1)*n_setv+1))
                   if(lsmallgrid) then
                      call interpolate(tmpreal2,tmpreal2,-1)
                   endif
                endif


                do iw=(iiw-1)*n_set+1,min(iiw*n_set,numw_prod)

                   if(.not.lsmallgrid) then

                      tmpspacec(:)=(0.d0,0.d0)
                      tmpspacec(nl(1:ngm))=tmpspacei(1:ngm,iw-(iiw-1)*n_set)
                      tmpspacec(nlm(1:ngm))=CONJG(tmpspacei(1:ngm,iw-(iiw-1)*n_set))

                      CALL invfft ('Dense', tmpspacec, dfftp)

                      tmpspacei1(:)=(0.d0,0.d0)
                      tmpspacec(:)=dble(tmpspacec(:)*tmpreal1(:))
                      CALL fwfft ('Dense', tmpspacec, dfftp)
                      tmpspacei1(1:ngm)=tmpspacec(nl(1:ngm))

                      if(  hw < min(hhw*n_setv,nbnd_v) ) then
                         tmpspacec(:)=(0.d0,0.d0)
                         tmpspacec(nl(1:ngm))=tmpspacei(1:ngm,iw-(iiw-1)*n_set)
                         tmpspacec(nlm(1:ngm))=CONJG(tmpspacei(1:ngm,iw-(iiw-1)*n_set))
                         CALL invfft ('Dense', tmpspacec, dfftp)
                         tmpspacec(:)=dble(tmpspacec(:)*tmpreal2(:))
                         CALL fwfft ('Dense', tmpspacec, dfftp)
                         tmpspacei2(1:ngm)=tmpspacec(nl(1:ngm))
                      endif

                   else
                      tmpspacec(:)=(0.d0,0.d0)
                      tmpspacec(nls(igk0(1:npw0)))=tmpspacei(1:npw0,iw-(iiw-1)*n_set)
                      CALL invfft ('Wave', tmpspacec, dffts)
                      tmpspacec(1:dffts%nnr)=dble(tmpspacec(1:dffts%nnr)*tmpreal1(1:dffts%nnr))
                      CALL fwfft ('Wave', tmpspacec, dffts)
                      tmpspacei1(1:npw0)=tmpspacec(nls(igk0(1:npw0)))
                      if(  hw < min(hhw*n_setv,nbnd_v) ) then
                          tmpspacec(:)=(0.d0,0.d0)
                          tmpspacec(nls(igk0(1:npw0)))=tmpspacei(1:npw0,iw-(iiw-1)*n_set)
                          CALL invfft ('Wave', tmpspacec, dffts)
                          tmpspacec(1:dffts%nnr)=dble(tmpspacec(1:dffts%nnr)*tmpreal2(1:dffts%nnr))
                          CALL fwfft ('Wave', tmpspacec, dffts)
                          tmpspacei2(1:npw0)=tmpspacec(nls(igk0(1:npw0)))
                       endif
                   endif

                   if(iiw==jjw) then
                      jw_begin=iw
                   else
                      jw_begin=(jjw-1)*n_set+1
                   endif
                   do jw=jw_begin,min(jjw*n_set,numw_prod)
                      if(ok_table(iw-(iiw-1)*n_set,jw-(jjw-1)*n_set) /= 0 ) then
                         wpwp_psi=0.d0
                         if(.not.lsmallgrid) then
                            do ig=1,ngm
                               wpwp_psi=wpwp_psi+ &
                                    & 2.d0*dble(conjg(tmpspacei1(ig))*tmpspacej(ig,jw-(jjw-1)*n_set))
                            enddo
                            if(gstart == 2) then
                               wpwp_psi=wpwp_psi-&
                                    &dble(conjg(tmpspacei1(1))*tmpspacej(1,jw-(jjw-1)*n_set))
                            endif
                         else
                            do ig=1,npw0
                               wpwp_psi=wpwp_psi+ &
                                    & 2.d0*dble(conjg(tmpspacei(ig,iw-(iiw-1)*n_set))*tmpspacej(ig,jw-(jjw-1)*n_set))
                            enddo
                            if(gstart == 2) then
                               wpwp_psi=wpwp_psi-&
                                    &dble(conjg(tmpspacei(1,iw-(iiw-1)*n_set))*tmpspacej(1,jw-(jjw-1)*n_set))
                            endif
                         endif
                         call mp_sum(wpwp_psi)
                         if(ionode) write(iunterm) ok_table(iw-(iiw-1)*n_set,jw-(jjw-1)*n_set),hw,wpwp_psi
                         !SCRIVI wpwp_psi su file
                         if(  hw < min(hhw*n_setv,nbnd_v) ) then
                            wpwp_psi=0.d0
                            if(.not.lsmallgrid) then
                               do ig=1,ngm
                                  wpwp_psi=wpwp_psi+ &
                                       & 2.d0*dble(conjg(tmpspacei2(ig))*tmpspacej(ig,jw-(jjw-1)*n_set))
                               enddo
                               if(gstart == 2) then
                                  wpwp_psi=wpwp_psi-&
                                       &dble(conjg(tmpspacei2(1))*tmpspacej(1,jw-(jjw-1)*n_set))
                               endif
                            else
                               do ig=1,npw0
                                  wpwp_psi=wpwp_psi+ &
                                       & 2.d0*dble(conjg(tmpspacei2(ig))*tmpspacej(ig,jw-(jjw-1)*n_set))
                               enddo
                               if(gstart == 2) then
                                  wpwp_psi=wpwp_psi-&
                                       &dble(conjg(tmpspacei2(1))*tmpspacej(1,jw-(jjw-1)*n_set))
                               endif
                            endif
                            call mp_sum(wpwp_psi)
                            if(ionode) write(iunterm) ok_table(iw-(iiw-1)*n_set,jw-(jjw-1)*n_set),hw+1,wpwp_psi
                         endif
                      endif
                   enddo
                enddo
             enddo
             deallocate(tmpspacej)
          enddo
          deallocate(tmpspacei)
       enddo
    enddo
    close(iungprod)
    if(ionode) then
       close(iunterm)
       close(iunterm2)
       allocate(index_terms(2,numw_prodprod))
       open( unit= iunterm2, file=trim(prefix)//'.wpwp_psi_index', status='old',action='readwrite',form='unformatted')
       rewind(iunterm2)
       read(iunterm2) iw!
       read(iunterm2) iw!nbnd_v
       read(iunterm2) iw!numw_prodprod
       do ig=1, numw_prodprod
          read(iunterm2) index_terms(1,ig), index_terms(2,ig)
       enddo
       rewind(iunterm2)
       write(iunterm2) numw_prod
       write(iunterm2) nbnd_v
       write(iunterm2) numw_prodprod
       do ig=1, numw_prodprod
          write(iunterm2) index_terms(1,ig), index_terms(2,ig)
       enddo
       write(stdout,*) 'TOTAL NUMBER OF PRODUCTS OF WANNIER PRODUCTS :', numw_prodprod
       close(iunterm2)
       deallocate(index_terms)
    endif
    deallocate(tmpspacev)
    deallocate(tmpreal1, tmpreals1)
    deallocate(tmpreal2, tmpreals2)
    deallocate(tmpspacec)
    deallocate(tmpspacei1,tmpspacei2)
    deallocate(tmpspaced)
    return
  end subroutine wannier_valence_terms_distance
