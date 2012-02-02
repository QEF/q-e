! FOR GWW
!
! Author: P. Umari
!
!-----------------------------------------------------------------------
subroutine matrix_wannier_gamma_big( matsincos, ispin, n_set, itask )
  !-----------------------------------------------------------------------
  !
  !this subroutine  calculates the terms <Psi_i|exp(iGX)|Psi_j>
  !in real space for gamma only case

! #ifdef __GWW

  USE kinds,                ONLY : DP
  USE cell_base,            ONLY : at, alat, tpiba, omega, tpiba2
  USE constants,            ONLY : e2, pi, tpi, fpi
  USE uspp,                 ONLY : okvan, nkb
  USE io_files,             ONLY : find_free_unit, diropn
  USE io_global,            ONLY : stdout
  USE realus,               ONLY : qsave, box,maxbox
  USE wannier_gw,           ONLY : becp_gw, expgsave, becp_gw_c, maxiter2,num_nbnd_first,num_nbndv,nbnd_normal
  USE ions_base,            ONLY : nat, ntyp => nsp, ityp
  USE uspp_param,           ONLY : upf, lmaxq, nh, nhm
  USE lsda_mod,             ONLY : nspin
  USE mp_global,            ONLY : intra_image_comm, me_pool
  USE fft_base,             ONLY : dfftp, dffts
  USE mp,                   ONLY : mp_bcast, mp_barrier, mp_sum

 implicit none
  !
  INTEGER, INTENT(in) :: ispin!spin polarization considred
!  COMPLEX(dp), INTENT(out) :: matp(nbnd_normal,nbnd_normal,3)
  REAL(dp), INTENT(out) :: matsincos(nbnd_normal,nbnd_normal,6)
  INTEGER, INTENT(in)  :: n_set  !defines the number of states
  INTEGER, INTENT(in)  :: itask !if ==1 consider subspace {C'}

  INTEGER ::  iiw,jjw, jw_begin
  INTEGER :: iw,jw,ir,ix,iy,iz,nn,ii
  REAL(kind=DP), ALLOCATABLE :: tmprealis(:,:),tmprealjs(:,:), tmpreal(:)
  COMPLEX(kind=DP), ALLOCATABLE ::  tmpexp(:), tmpexp2(:,:)
  INTEGER ::  iunwfcreal2
  COMPLEX(kind=DP) :: sca,ee,sca1
  REAL(kind=DP) :: dsgn
  LOGICAL :: exst
  INTEGER :: iqq
  INTEGER :: nt, na, ih, jh, np, mbia, irb, iqs, nhnt, ia
  INTEGER :: ikb, jkb, ijkb0, is
  INTEGER :: isgn,mdir
  INTEGER :: nr3s_start, nr3s_end
  INTEGER :: nr3_start, nr3_end
  INTEGER :: nbnd_eff

  write(stdout,*) 'MATRIX BIG1'
  call flush_unit(stdout)

  iunwfcreal2=find_free_unit()
  CALL diropn( iunwfcreal2, 'real_whole', dffts%nnr, exst )


  allocate(tmprealis(dffts%nnr,n_set),tmprealjs(dffts%nnr,n_set), tmpreal(dffts%nnr))
  allocate(tmpexp2(dffts%nnr,6))

!set up exponential grid

  tmpexp2(:,:)=(0.d0,0.d0)

#ifndef __MPI
  iqq=0
  do ix=1,dffts%nr1
     do iy=1,dffts%nr2
        do iz=1,dffts%nr3
           iqq=(iz-1)*(dffts%nr1x*dffts%nr2x)+(iy-1)*dffts%nr1x+ix
           tmpexp2(iqq,1)= exp(cmplx(0.d0, 1.d0)*tpi*real(ix-1)/real(dffts%nr1))
           tmpexp2(iqq,2)= exp(cmplx(0.d0, 1.d0)*tpi*real(iy-1)/real(dffts%nr2))
           tmpexp2(iqq,3)= exp(cmplx(0.d0, 1.d0)*tpi*real(iz-1)/real(dffts%nr3))
           tmpexp2(iqq,4)= exp(cmplx(0.d0,-1.d0)*tpi*real(ix-1)/real(dffts%nr1))
           tmpexp2(iqq,5)= exp(cmplx(0.d0,-1.d0)*tpi*real(iy-1)/real(dffts%nr2))
           tmpexp2(iqq,6)= exp(cmplx(0.d0,-1.d0)*tpi*real(iz-1)/real(dffts%nr3))
        enddo
     enddo
  enddo


#else
  write(stdout,*) 'NRS' , dffts%nr1 ,dffts%nr2 ,dffts%nr3
  write(stdout,*) 'NRXS', dffts%nr1x,dffts%nr2x,dffts%nr3x
  nr3s_start=0
  nr3s_end =0
  do ii=1,me_pool + 1
     nr3s_start=nr3s_end+1
     nr3s_end=nr3s_end+dffts%npp(ii)
  end do
  tmpexp2(:,:)=(0.d0,0.d0)
  do iz=1,dffts%npp(me_pool+1)
     do iy=1,dffts%nr2
        do ix=1,dffts%nr1
           iqq=(iz-1)*(dffts%nr1x*dffts%nr2x)+(iy-1)*dffts%nr1x+ix
           tmpexp2(iqq,1) = exp(cmplx(0.d0,1.d0)*tpi*real(ix-1)/real(dffts%nr1))
           tmpexp2(iqq,2) = exp(cmplx(0.d0,1.d0)*tpi*real(iy-1)/real(dffts%nr2))
           tmpexp2(iqq,3) = exp(cmplx(0.d0,1.d0)*tpi*real(iz+nr3s_start-1-1)/real(dffts%nr3))
           tmpexp2(iqq,4) = exp(cmplx(0.d0,-1.d0)*tpi*real(ix-1)/real(dffts%nr1))
           tmpexp2(iqq,5) = exp(cmplx(0.d0,-1.d0)*tpi*real(iy-1)/real(dffts%nr2))
           tmpexp2(iqq,6) = exp(cmplx(0.d0,-1.d0)*tpi*real(iz+nr3s_start-1-1)/real(dffts%nr3))
        enddo
     enddo
  enddo


#endif

  write(stdout,*) 'Calculate grid'


  if(maxiter2 >= 1 .or. num_nbnd_first==0) then
     nbnd_eff=nbnd_normal
  else
     nbnd_eff=num_nbndv+num_nbnd_first
  endif

  write(stdout,*) 'MATRIX BIG2'
  call flush_unit(stdout)

  do iiw=1,nbnd_eff/n_set+1
     write(stdout,*) 'MATRIX IIW',iiw
     call flush_unit(stdout)

     do iw=(iiw-1)*n_set+1,min(iiw*n_set,nbnd_eff)
!read from disk wfc on coarse grid
        CALL davcio( tmprealis(:,iw-(iiw-1)*n_set),dffts%nnr,iunwfcreal2,iw,-1)
     enddo
!read in iw wfcs
     do jjw=iiw,nbnd_eff/n_set+1
        write(stdout,*) 'MATRIX JJW',jjw
        call flush_unit(stdout)

        do jw=(jjw-1)*n_set+1,min(jjw*n_set,nbnd_eff)
           CALL davcio( tmprealjs(:,jw-(jjw-1)*n_set),dffts%nnr,iunwfcreal2,jw,-1)
        enddo
        !do product

        do iw=(iiw-1)*n_set+1,min(iiw*n_set,nbnd_eff)
           if(iiw==jjw) then
              jw_begin=iw
           else
              jw_begin=(jjw-1)*n_set+1
           endif
           do jw=jw_begin,min(jjw*n_set,nbnd_eff)

              tmpreal(:)=tmprealis(:,iw-(iiw-1)*n_set)*tmprealjs(:,jw-(jjw-1)*n_set)

!put on fine grid


!add us part

!calculate matrix element
              do mdir=1,3
                 sca=0.d0
                 do ir=1,dffts%nnr
                    sca=sca+tmpreal(ir)*tmpexp2(ir,mdir)
                 enddo
                 sca=sca/dble(dffts%nr1*dffts%nr2*dffts%nr3)
                 call mp_barrier
                 call mp_sum(sca)
                 matsincos(iw,jw,mdir)=dble(sca)
                 matsincos(jw,iw,mdir)=dble(sca)
                 matsincos(iw,jw,mdir+3)=dimag(sca)
                 matsincos(jw,iw,mdir+3)=dimag(sca)
                 !matp(iw,jw,mdir)=sca
                 !matp(jw,iw,mdir)=sca
              enddo


           enddo
        enddo
     enddo
  enddo


  deallocate(tmprealis,tmprealjs)
  deallocate(tmpexp2)

  write(stdout,*) 'Calculate US'
  call flush_unit(stdout)
  if(okvan) then
    allocate(tmpexp(dfftp%nnr))
    allocate(expgsave(maxval(nh),maxval(nh),nat,3))
    expgsave(:,:,:,:)=0.d0
   do mdir=1,3

#ifndef __MPI
      if(mdir==1) then
         do ix=1,dfftp%nr1
            ee=exp(cmplx(0.d0,1.d0)*tpi*real(ix)/real(dfftp%nr1))
            do iy=1,dfftp%nr2
              do  iz=1,dfftp%nr3
                 nn=(iz-1)*dfftp%nr1x*dfftp%nr2x+(iy-1)*dfftp%nr1x+ix
                 tmpexp(nn)=ee
              enddo
           enddo
         enddo
      else if(mdir==2) then
         do iy=1,dfftp%nr2
            ee=exp(cmplx(0.d0,1.d0)*tpi*real(iy)/real(dfftp%nr2))
            do ix=1,dfftp%nr1
              do  iz=1,dfftp%nr3
                 nn=(iz-1)*dfftp%nr1x*dfftp%nr2x+(iy-1)*dfftp%nr1x+ix
                 tmpexp(nn)=ee
              enddo
            enddo
         enddo
      else if(mdir==3) then
         do iz=1,dfftp%nr3
         ee=exp(cmplx(0.d0,1.d0)*tpi*real(iz)/real(dfftp%nr3))
            do ix=1,dfftp%nr1
              do  iy=1,dfftp%nr2
                 nn=(iz-1)*dfftp%nr1x*dfftp%nr2x+(iy-1)*dfftp%nr1x+ix
                 tmpexp(nn)=ee
              enddo
            enddo
         enddo
      endif

#else
      nr3_start=0
      nr3_end =0
      do ii=1,me_pool + 1
         nr3_start=nr3_end+1
         nr3_end=nr3_end+dfftp%npp(ii)
      end do

      do iz=1,dfftp%npp(me_pool+1)
         do iy=1,dfftp%nr2
            do ix=1,dfftp%nr1

               nn=(iz-1)*dfftp%nr1x*dfftp%nr2x+(iy-1)*dfftp%nr1x+ix
               if(mdir==1) then
                  tmpexp(nn)= exp(cmplx(0.d0,1.d0)*tpi*real(ix-1)/real(dfftp%nr1))
               elseif(mdir==2) then
                  tmpexp(nn)= exp(cmplx(0.d0,1.d0)*tpi*real(iy-1)/real(dfftp%nr2))
               else
                  tmpexp(nn)= exp(cmplx(0.d0,1.d0)*tpi*real(iz+nr3_start-1-1)/real(dfftp%nr3))
               endif
            enddo
         enddo
      enddo


#endif

     ijkb0 = 0
     DO np = 1, ntyp
        !
        iqs = 0
        !
        IF ( upf(np)%tvanp ) THEN
           !
           DO ia = 1, nat
              !
              mbia = maxbox(ia)
              !
              nt = ityp(ia)
              nhnt = nh(nt)
              !
              IF ( ityp(ia) /= np ) iqs=iqs+(nhnt+1)*nhnt*mbia/2
              IF ( ityp(ia) /= np ) CYCLE
              !
              DO ih = 1, nhnt
                 !
                 ikb = ijkb0 + ih
                 !
                 DO jh = ih, nhnt
                    !
                    jkb = ijkb0 + jh
                    !
                    expgsave(ih,jh,ia,mdir)=(0.d0,0.d0)
                    DO ir = 1, mbia
                       !
                       irb = box(ir,ia)
                       iqs = iqs + 1
                       !
                       expgsave(ih,jh,ia,mdir)=expgsave(ih,jh,ia,mdir)+qsave(iqs)*tmpexp(irb)
                       !
                    ENDDO
                 ENDDO
              ENDDO
              ijkb0 = ijkb0 + nhnt
              !
           ENDDO
           !
        ELSE
           !
           DO ia = 1, nat
              !
              IF ( ityp(ia) == np ) ijkb0 = ijkb0 + nh(np)
              !
           END DO
           !
        END IF
     ENDDO

     expgsave(:,:,:,mdir)=expgsave(:,:,:,mdir)*omega/dble(dfftp%nr1*dfftp%nr2*dfftp%nr3)

#ifdef __MPI
     !!!call reduce (2  *maxval(nh) *maxval(nh)* nat, expgsave(:,:,:,mdir))
     call mp_sum(expgsave(:,:,:,mdir))
#endif


      do iw=1,nbnd_eff
        do jw=iw,nbnd_eff

          do is=1, nspin
           ijkb0 = 0
           do  np = 1, ntyp
             if ( upf(np)%tvanp ) then
               do  na = 1, nat
                 if ( ityp(na) == np ) then
                  do ih = 1, nh(np)
                    ikb = ijkb0 + ih
                    do jh = 1, nh(np)
                      jkb = ijkb0 + jh
                      if(ih <= jh) then
                         if(itask /= 1) then
                            !matp(iw,jw,mdir)=matp(iw,jw,mdir)+expgsave(ih,jh,na,mdir) * becp_gw(ikb,iw)*becp_gw(jkb,jw)
                            matsincos(iw,jw,mdir)=matsincos(iw,jw,mdir)+&
                                 &dble(expgsave(ih,jh,na,mdir) * becp_gw(ikb,iw)*becp_gw(jkb,jw))
                             matsincos(iw,jw,mdir+3)=matsincos(iw,jw,mdir+3)+&
                                 &dimag(expgsave(ih,jh,na,mdir) * becp_gw(ikb,iw)*becp_gw(jkb,jw))
                         else
                            !matp(iw,jw,mdir)=matp(iw,jw,mdir)+expgsave(ih,jh,na,mdir) * becp_gw_c(ikb,iw)*becp_gw_c(jkb,jw)
                            matsincos(iw,jw,mdir)=matsincos(iw,jw,mdir)+&
                                 &dble(expgsave(ih,jh,na,mdir) * becp_gw_c(ikb,iw)*becp_gw_c(jkb,jw))
                            matsincos(iw,jw,mdir+3)=matsincos(iw,jw,mdir+3)+&
                                 &dimag(expgsave(ih,jh,na,mdir) * becp_gw_c(ikb,iw)*becp_gw_c(jkb,jw))
                        endif

                      else
                         if(itask /= 1) then
                            !matp(iw,jw,mdir)=matp(iw,jw,mdir)+expgsave(jh,ih,na,mdir)  * becp_gw(ikb,iw)*becp_gw(jkb,jw)
                            matsincos(iw,jw,mdir)=matsincos(iw,jw,mdir)+&
                                 &dble(expgsave(jh,ih,na,mdir)  * becp_gw(ikb,iw)*becp_gw(jkb,jw))
                             matsincos(iw,jw,mdir+3)=matsincos(iw,jw,mdir+3)+&
                                 &dimag(expgsave(jh,ih,na,mdir)  * becp_gw(ikb,iw)*becp_gw(jkb,jw))
                         else
                            !matp(iw,jw,mdir)=matp(iw,jw,mdir)+expgsave(jh,ih,na,mdir)  * becp_gw_c(ikb,iw)*becp_gw_c(jkb,jw)
                            matsincos(iw,jw,mdir)=matsincos(iw,jw,mdir)+&
                                 &dble(expgsave(jh,ih,na,mdir)  * becp_gw_c(ikb,iw)*becp_gw_c(jkb,jw))
                             matsincos(iw,jw,mdir+3)=matsincos(iw,jw,mdir+3)+&
                                 &dimag(expgsave(jh,ih,na,mdir)  * becp_gw_c(ikb,iw)*becp_gw_c(jkb,jw))
                         endif
                      endif
                    enddo
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
         enddo
!        matp(jw,iw,mdir)=matp(iw,jw,mdir)
         matsincos(jw,iw,mdir)=matsincos(iw,jw,mdir)
         matsincos(jw,iw,mdir+3)=matsincos(iw,jw,mdir+3)
    enddo
   enddo
  enddo

  deallocate(tmpexp)
  endif

  close(iunwfcreal2)
! #endif

  return

  end subroutine


