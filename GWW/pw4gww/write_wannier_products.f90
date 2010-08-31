! FOR GWW
!
! Author: P. Umari
!
   subroutine write_wannier_products( wp, ifile,num_wp_dimr,num_wp_r)
!this subroutine writes the data relative
!to the overlaps of orthonormalized products of wannier
!with products of wannier
!to be read by GWW code


  USE wannier_gw, ONLY : numw_prod, wannier_P, l_pmatrix
  USE io_files, ONLY : find_free_unit, prefix
  USE mp_global, ONLY : mpime, nproc
  USE io_global, ONLY : ionode, ionode_id
  USE mp, ONLY : mp_bcast
  USE kinds, ONLY : DP

  implicit none

  TYPE(wannier_P), INTENT(in) ::   wp(num_wp_dimr)!products of wanniers
  INTEGER, INTENT(in) :: ifile!0=file of products 1=file of terms w_i(r)w_j(r)v(r,r')w^P_red(r')
  INTEGER, INTENT(in) :: num_wp_dimr!dimension of wp
  INTEGER, INTENT(in) :: num_wp_r!periodicity of qp for parallel execution

  INTEGER :: iunu, iw
#ifdef __SCALAPACK
  INTEGER, EXTERNAL :: indxg2l,indxg2p
#endif
  INTEGER :: numij
  INTEGER, ALLOCATABLE :: ij(:,:)
  REAL(kind=DP), ALLOCATABLE :: o(:)
  INTEGER :: icrow,ilrow

  if(ionode) then
     iunu = find_free_unit()

     if(ifile==0) then
        open(unit=iunu,file=trim(prefix)//'.wp',status='unknown',form='unformatted')
     else
        open(unit=iunu,file=trim(prefix)//'.wp_v',status='unknown',form='unformatted')
     endif

     write(iunu) numw_prod
  endif

  if(.not.l_pmatrix) then
     if(ionode) then
        do iw=1,numw_prod
           write(iunu) wp(iw)%numij
           write(iunu) wp(iw)%ij(1,1:wp(iw)%numij)
           write(iunu) wp(iw)%ij(2,1:wp(iw)%numij)
           write(iunu) wp(iw)%o(1:wp(iw)%numij)
        enddo
     endif
  else
#ifdef __SCALAPACK
     do iw=1,numw_prod
        icrow=indxg2p(iw,num_wp_r,0,0,nproc)
        if(icrow==mpime) then
           ilrow=indxg2l(iw,num_wp_r,0,0,nproc)
           numij=wp(ilrow)%numij
        endif
        call mp_bcast(numij,icrow)
        allocate(ij(2,numij))
        allocate(o(numij))
        if(icrow==mpime) then
           ij(1:2,1:numij) =wp(ilrow)%ij(1:2,1:numij)
           o(1:numij)=wp(ilrow)%o(1:numij)
        endif
        call mp_bcast(ij,icrow)
        call mp_bcast(o,icrow)
        if(ionode) then
           write(iunu) numij
           write(iunu) ij(1,1:numij)
           write(iunu) ij(2,1:numij)
           write(iunu) o(1:numij)
        endif
        deallocate(ij,o)
     enddo
#endif
  endif

  if(ionode) close(iunu)

  return
  end subroutine
  !
  ! -----------------------------------------------
  !
  subroutine orthonormalize_producs_cutoff(numpw,numpw_max,numpw_offset,numpw_original,&
                       omat, numpwdimr,numpwdimc,cutoff, n_set, ecutoff, lcutoff)
!this subroutine calculate the inverse transformation matrix U
!\{w}^P_i = U_ij \tilde{w}^P_j
!then apply a cut off on the eigenvalues
!and write the orthonormalize wannier products on file
!Note:for calculations without scalapack omat_eff is now local to the first processor


     USE kinds, ONLY : DP
     USE io_global, ONLY : stdout, ionode, ionode_id
     USE io_files, ONLY : find_free_unit, prefix, diropn
     USE gvect,  ONLY : ngm, gg
     USE cell_base, ONLY: tpiba2
     USE mp,        ONLY : mp_bcast, mp_sum,mp_barrier
     USE wannier_gw, ONLY : max_ngm,p_mpime,p_nproc, npcol, nprow,icontxt,myrow,mycol,l_pmatrix
     USE mp_global, ONLY : intra_pool_comm

     implicit none

     INTEGER, INTENT(inout) :: numpw!dimension of matrices in  MODIFIED IN OUTPUT
     INTEGER, INTENT(in) :: numpw_max!maximum number on which orthogonalizing
     INTEGER, INTENT(in) :: numpw_offset!offset for first product to be considered, starting from 0
!it is orthogonalizing from numpw_offset+1 and numpw_max
     INTEGER, INTENT(out) :: numpw_original!dimension of original matrices
     REAL(kind=DP) :: omat(numpwdimr,numpwdimc)!overlap matrix ON OUTPUT contains the terms <\tilde{w}_j|w_i>
     INTEGER,INTENT(in) :: numpwdimr,numpwdimc!dimensions of omat
     REAL(kind=DP) :: cutoff!cutoff for eigenvalues
     INTEGER :: n_set !length of buffer for products of wannier
     REAL(kind=DP), INTENT(in) :: ecutoff!cutoff in Rydberg for g sum
     LOGICAL, INTENT(in) :: lcutoff !if true uses cutoff on G defined by ecutoff POTENTIALLY DANGEROUS



     REAL(kind=DP), ALLOCATABLE :: eigen(:),work(:)
     INTEGER :: lwork,info,liwork
     INTEGER, ALLOCATABLE :: iwork(:)

     INTEGER :: i,j,k, ig, ip, ii
     REAL(kind=DP) :: fact
     INTEGER :: iunu, iungprod,iungprod_on
     COMPLEX(kind=DP), ALLOCATABLE :: tmpspacei(:,:)
     COMPLEX(kind=DP), ALLOCATABLE :: tmpspacej(:,:)
     INTEGER :: iw,jw,iiw,jjw
     INTEGER, ALLOCATABLE :: on_cor(:)
     INTEGER :: numpw_on
     REAL(kind=DP), ALLOCATABLE :: tmp_str(:)
     LOGICAL :: exst
     INTEGER :: ngm_max
     REAL(kind=DP), ALLOCATABLE :: omat_eff(:,:), omat_tmp(:,:)
     INTEGER :: numpw_eff


     INTEGER :: numpw_eff_r,numpw_eff_c
     INTEGER :: numpw_r,numpw_c,numpw_eff_dimr,numpw_eff_dimc,numpw_dimr,numpw_dimc
     REAL(kind=DP), ALLOCATABLE :: omat_eff_l(:,:), omat_l(:,:)
     INTEGER :: desc_a(9),desc_b(9),desc_c(9)
     INTEGER :: numpw_on_r,numpw_on_c,numpw_on_dimr,numpw_on_dimc
     INTEGER, ALLOCATABLE :: ifail(:),iclustr(:)
     REAL(kind=DP), ALLOCATABLE :: gap(:)
     INTEGER :: m,nz,icrow,iccol,iproc,ilrow,ilcol
#ifdef __SCALAPACK
     INTEGER, EXTERNAL :: indxg2p,indxg2l
#endif
     REAL(kind=DP) :: rsca

     REAL(kind=DP), ALLOCATABLE :: sndbuf(:,:)
     INTEGER, ALLOCATABLE :: ileng(:), ipos(:),ibuf(:,:)
     REAL(kind=DP), ALLOCATABLE :: omat_eff_bcast(:) !part of omat_eff to be broadcast form 1st processor


     call mp_barrier!ATTENZIONE
     write(stdout,*) 'ROUTINE orthonormalize_producs_cutoff',numpw,numpw_max,numpw_offset,numpw_original,&
          &numpwdimr,numpwdimc,cutoff, n_set, ecutoff, lcutoff
     call flush_unit(stdout)

     numpw_eff=numpw_max-numpw_offset

     !determine ngm_max
     if(lcutoff) then
        ngm_max=0
        do ig=1,ngm
           if(gg(ig)*tpiba2 >= ecutoff) exit
           ngm_max=ngm_max+1
        enddo
     else
        ngm_max=ngm
     endif

     write(stdout,*) 'NGM MAX:', ngm_max, ngm
     call flush_unit(stdout)


     if(.not.l_pmatrix) then
        allocate(eigen(numpw_eff))
        allocate(omat_eff_bcast(numpw_eff))
        if(ionode) then
           allocate(omat_eff(numpw_eff,numpw_eff))
           omat_eff(1:numpw_eff,1:numpw_eff)=&
                &omat(numpw_offset+1:numpw_max,numpw_offset+1:numpw_max)
           allocate(work(1))
           call  DSYEV( 'V', 'U', numpw_eff, omat_eff, numpw_eff, eigen, work, -1, info )
           lwork=work(1)
           deallocate(work)
           allocate(work(lwork))
           call  DSYEV( 'V', 'U', numpw_eff, omat_eff, numpw_eff, eigen, work, lwork, info )
           deallocate(work)
           if(info/=0) then
              write(stdout,*) 'ROUTINE find_transform, INFO:', info
              stop
           endif
        endif
!        do iw=1,numpw_eff
!           call mp_barrier
!           call mp_bcast(omat_eff(:,iw), ionode_id)
!        enddo
        call mp_bcast(eigen(:), ionode_id)
     else
#ifdef __SCALAPACK

        write(stdout,*) 'ORTHO1'
        call flush_unit(stdout)

!set up omat_eff
        numpw_eff_r=ceiling(real(numpw_eff)/real(max(nprow,npcol)))
        numpw_eff_c=ceiling(real(numpw_eff)/real(max(nprow,npcol)))
    !    numpw_r=ceiling(real(numpw)/real(nprow))
    !    numpw_c=ceiling(real(numpw)/real(npcol))
        numpw_r=ceiling(real(numpw)/real(max(npcol,nprow)))
        numpw_c=ceiling(real(numpw)/real(max(nprow,npcol)))
!calculate dimesion of local omat_tmp and allocate
        numpw_eff_dimr=ceiling (real(numpw_eff)/real(numpw_eff_r*nprow))*numpw_eff_r
        numpw_eff_dimc=ceiling (real(numpw_eff)/real(numpw_eff_c*npcol))*numpw_eff_c
        numpw_dimr=ceiling (real(numpw)/real(numpw_r*nprow))*numpw_r
        numpw_dimc=ceiling (real(numpw)/real(numpw_c*npcol))*numpw_c
        if((myrow+1) <= nprow .and. (mycol+1) <= npcol) then
           allocate(omat_eff_l(numpw_eff_dimr,numpw_eff_dimc))
           allocate(omat_eff(numpw_eff_dimr,numpw_eff_dimc))
           allocate(eigen(numpw_eff))
           allocate(gap(numpw_eff),ifail(numpw_eff),iclustr(numpw_eff))

           write(stdout,*) 'ORTHO2'
           call flush_unit(stdout)


           allocate(sndbuf(numpw_eff,nprow*npcol),ibuf(numpw_eff,nprow*npcol))
           allocate(ipos(nprow*npcol))

            do jw=1,numpw_eff
               write(stdout,*) 'ORTHO2.5',jw
               call flush_unit(stdout)

               ipos(:)=0

               do iw=1,numpw_eff
                  icrow = indxg2p(jw+numpw_offset,numpw_r,0,0,nprow)
                  iccol = indxg2p(iw+numpw_offset,numpw_c,0,0,npcol)
                  iproc=icrow*npcol+iccol
                  ipos(iproc+1)=ipos(iproc+1)+1
                  ibuf(ipos(iproc+1),iproc+1)=iw
                  if(myrow==icrow .and. mycol==iccol) then
                     ilrow=indxg2l(jw+numpw_offset,numpw_r,0,0,nprow)
                     ilcol=indxg2l(iw+numpw_offset,numpw_c,0,0,npcol)
                     sndbuf(ipos(iproc+1),iproc+1)=omat(ilrow,ilcol)
                  endif
               enddo
               do ip=1,nprow*npcol
                  if(ipos(ip)>0) then
                     !call mp_barrier
                     call mp_bcast(sndbuf(1:ipos(ip),ip),ip-1)
                     do ii=1,ipos(ip)
                        iw=ibuf(ii,ip)
                        icrow = indxg2p(jw,numpw_eff_r,0,0,nprow)
                        iccol = indxg2p(iw,numpw_eff_c,0,0,npcol)
                         if(myrow==icrow .and. mycol==iccol) then
                            ilrow=indxg2l(jw,numpw_eff_r,0,0,nprow)
                            ilcol=indxg2l(iw,numpw_eff_c,0,0,npcol)
                            omat_eff_l(ilrow,ilcol)=sndbuf(ii,ip)
                         endif
                      enddo
                   endif
                enddo



!                  icrow = indxg2p(jw+numpw_offset,numpw_r,0,0,nprow)
!                  iccol = indxg2p(iw+numpw_offset,numpw_c,0,0,npcol)
!                  iproc=icrow*npcol+iccol
!                  if(myrow==icrow .and. mycol==iccol) then
!                     ilrow=indxg2l(jw+numpw_offset,numpw_r,0,0,nprow)
!                     ilcol=indxg2l(iw+numpw_offset,numpw_c,0,0,npcol)
!                     rsca=omat(ilrow,ilcol)
!                  endif
!                  call mp_bcast(rsca,iproc)
!                  icrow = indxg2p(jw,numpw_eff_r,0,0,nprow)
!                  iccol = indxg2p(iw,numpw_eff_c,0,0,npcol)
!                  if(myrow==icrow .and. mycol==iccol) then
!                     ilrow=indxg2l(jw,numpw_eff_r,0,0,nprow)
!                     ilcol=indxg2l(iw,numpw_eff_c,0,0,npcol)
!                      omat_eff_l(ilrow,ilcol)=rsca
!                  endif
             enddo

            deallocate(sndbuf,ibuf)
            deallocate(ipos)


!        omat_eff(1:numpw_eff,1:numpw_eff)=omat(numpw_offset+1:numpw_max,numpw_offset+1:numpw_max)
! A = omat_eff_l
             desc_a(1)=1
             desc_a(2)=icontxt
             desc_a(3)=numpw_eff
             desc_a(4)=numpw_eff
             desc_a(5)=numpw_eff_r
             desc_a(6)=numpw_eff_c
             desc_a(7)=0
             desc_a(8)=0
             desc_a(9)=numpw_eff_dimr
!B= omat_eff
             desc_b(1)=1
             desc_b(2)=icontxt
             desc_b(3)=numpw_eff
             desc_b(4)=numpw_eff
             desc_b(5)=numpw_eff_r
             desc_b(6)=numpw_eff_c
             desc_b(7)=0
             desc_b(8)=0
             desc_b(9)=numpw_eff_dimr

             allocate(work(1),iwork(1))

             write(stdout,*) 'ORTHO3'
             call flush_unit(stdout)




!             call pdsyevd('V','U',numpw_eff,omat_eff_l,1,1,desc_a,eigen,omat_eff,1,1,desc_b,work,-1,iwork,-1,info)
             call pdsyevx('V','A','U',numpw_eff,omat_eff_l,1,1,desc_a,0.d0,0.d0,0,0,0.d0,&
 & m,nz,eigen,-1.d0,omat_eff,1,1,desc_b,work,-1,iwork,-1,ifail,iclustr,gap,info)
             lwork=work(1)
             liwork=iwork(1)
             deallocate(work,iwork)
             write(stdout,*) 'LWORK LIWORK:',lwork,liwork
             call flush_unit(stdout)
             allocate(work(lwork),iwork(liwork))
!             call pdsyevd('V','U',numpw_eff,omat_eff_l,1,1,desc_a,eigen,omat_eff,1,1,desc_b,work,lwork,iwork,liwork,info)
               call pdsyevx('V','A','U',numpw_eff,omat_eff_l,1,1,desc_a,0.d0,0.d0,0,0,0.d0,&
       & m,nz,eigen,0.d0,omat_eff,1,1,desc_b,work,lwork,iwork,liwork,ifail,iclustr,gap,info)
             if(info/=0) then
                write(stdout,*) 'ROUTINE PDSYEV INFO:', info
                stop
             endif


             deallocate(work,iwork)
             deallocate(gap,iclustr,ifail)
             deallocate(omat_eff_l)

!  call  DSYEV( 'V', 'U', numpw_eff, omat_eff, numpw_eff, eigen, work, -1, info )

          endif
#endif
      endif


      call mp_barrier
      write(stdout,*) 'out of diagonalization'
      call flush_unit(stdout)
      call mp_bcast(eigen(:),ionode_id)

     if(ionode) write(stdout,*) 'OMAT EFF', omat_eff(1,1)!ATTENZIONE

     do iw=1,numpw_eff
        write(stdout,*) 'EIGEN:',iw, eigen(iw)
     enddo



     iungprod = find_free_unit()
     CALL diropn( iungprod, 'wiwjwfc', max_ngm*2, exst )
     iungprod_on = find_free_unit()
     CALL diropn( iungprod_on, 'wiwjwfc_on', max_ngm*2, exst )



!allocate file
     allocate(tmpspacei(max_ngm,n_set),tmpspacej(max_ngm,n_set))
     allocate(on_cor(numpw_eff))
     on_cor(:)=0

     numpw_on=0
!external loop on orthonormal products
     do iiw=1,ceiling(real(numpw_eff)/real(n_set))
!if required read them from disk or set to zero
        write(stdout,*) 'IIW', iiw
        call flush_unit(stdout)
        tmpspacei(:,:)=0.d0
        do jjw=1,ceiling(real(numpw_eff)/real(n_set))
            write(stdout,*) 'JJW', jjw
            call flush_unit(stdout)
            do jw=(jjw-1)*n_set+1,min(jjw*n_set,numpw_eff)
               CALL davcio(tmpspacej(:,jw-(jjw-1)*n_set), max_ngm*2,iungprod,jw+numpw_offset,-1)
            enddo
            do iw=(iiw-1)*n_set+1,min(iiw*n_set,numpw_eff)
               if(eigen(iw) >= cutoff) then
                  if(on_cor(iw) == 0) then
                     numpw_on=numpw_on+1
                     on_cor(iw)=numpw_on
                  endif

                  if(.not.l_pmatrix) then
                     if(ionode) omat_eff_bcast(:)=omat_eff(:,iw)
                     call mp_barrier
                     call mp_bcast(omat_eff_bcast, ionode_id)
                     do jw=(jjw-1)*n_set+1,min(jjw*n_set,numpw_eff)
!                        tmpspacei(1:ngm_max,iw-(iiw-1)*n_set)=tmpspacei(1:ngm_max,iw-(iiw-1)*n_set)+&
!                             &(1.d0/dsqrt(eigen(iw)))*omat_eff(jw,iw)*tmpspacej(1:ngm_max,jw-(jjw-1)*n_set)
                        tmpspacei(1:ngm_max,iw-(iiw-1)*n_set)=tmpspacei(1:ngm_max,iw-(iiw-1)*n_set)+&
                             &(1.d0/dsqrt(eigen(iw)))*omat_eff_bcast(jw)*tmpspacej(1:ngm_max,jw-(jjw-1)*n_set)


                     enddo
                  else
#ifdef __SCALAPACK
                     allocate(sndbuf(numpw_eff,nprow*npcol),ibuf(numpw_eff,nprow*npcol))
                     allocate(ipos(nprow*npcol))
                     ipos(:)=0

                     do jw=(jjw-1)*n_set+1,min(jjw*n_set,numpw_eff)

                        icrow = indxg2p(jw,numpw_eff_r,0,0,nprow)
                        iccol = indxg2p(iw,numpw_eff_c,0,0,npcol)
                        iproc=icrow*npcol+iccol
                        ipos(iproc+1)=ipos(iproc+1)+1
                        ibuf(ipos(iproc+1),iproc+1)=jw
                        if(myrow==icrow .and. mycol==iccol) then
                           ilrow=indxg2l(jw,numpw_eff_r,0,0,nprow)
                           ilcol=indxg2l(iw,numpw_eff_c,0,0,npcol)
                           sndbuf(ipos(iproc+1),iproc+1)=omat_eff(ilrow,ilcol)
                        endif

                     enddo
                     do ip=1,nprow*npcol
                        if(ipos(ip)>0) then
                           !call mp_barrier
                           call mp_bcast(sndbuf(1:ipos(ip),ip),ip-1)
                           do ii=1,ipos(ip)
                              jw=ibuf(ii,ip)
                              tmpspacei(1:ngm_max,iw-(iiw-1)*n_set)=tmpspacei(1:ngm_max,iw-(iiw-1)*n_set)+&
                                   &(1.d0/dsqrt(eigen(iw)))*sndbuf(ii,ip)*tmpspacej(1:ngm_max,jw-(jjw-1)*n_set)
                           enddo
                        endif
                     enddo

                     deallocate(sndbuf,ibuf,ipos)
#endif
                  endif
!                  do jw=(jjw-1)*n_set+1,min(jjw*n_set,numpw_eff)
!                     if(.not.l_pmatrix) then
!                        tmpspacei(1:ngm_max,iw-(iiw-1)*n_set)=tmpspacei(1:ngm_max,iw-(iiw-1)*n_set)+(1.d0/dsqrt(eigen(iw)))*omat_ef
!                             &tmpspacej(1:ngm_max,jw-(jjw-1)*n_set)
!                     else
!#ifdef __SCALAPACK
!determine if jw and iw stays in the processor
!USARE LAPACK ROUTINES
!                        icrow = indxg2p(jw,numpw_eff_r,0,0,nprow)
!                        iccol = indxg2p(iw,numpw_eff_c,0,0,npcol)
!                        iproc = icrow*npcol+iccol
!                        if(myrow==icrow .and. mycol==iccol) then
!                           ilrow=indxg2l(jw,numpw_eff_r,0,0,nprow)
!                           ilcol=indxg2l(iw,numpw_eff_c,0,0,npcol)
!                           rsca=omat_eff(ilrow,ilcol)
!                        endif
!                        call mp_bcast(rsca,iproc)


!                          tmpspacei(1:ngm_max,iw-(iiw-1)*n_set)=tmpspacei(1:ngm_max,iw-(iiw-1)*n_set)+&
! &(1.d0/dsqrt(eigen(iw)))*rsca*tmpspacej(1:ngm_max,jw-(jjw-1)*n_set)
!#endif
!                       endif
!                    enddo
               endif
            enddo
         enddo
         do iw=(iiw-1)*n_set+1,min(iiw*n_set,numpw_eff)
            if(on_cor(iw) /= 0) then
               CALL davcio(tmpspacei(:,iw-(iiw-1)*n_set), max_ngm*2,iungprod_on,on_cor(iw),1)
            endif
         enddo
      enddo




     close(iungprod)
     close(iungprod_on)


!for compatibility write orthonorm file


     if(ionode) then
        allocate(tmp_str(numpw_on))
        iunu = find_free_unit()

        open(unit=iunu,file=trim(prefix)//'.orthonorm',status='unknown',form='unformatted')
        write(iunu) numpw_on

        do  i=1,numpw_on
           tmp_str(:)=0.d0
           tmp_str(i)=1.d0
           write(iunu) tmp_str
        enddo

        close(iunu)
        deallocate(tmp_str)
     endif

     call mp_barrier
     write(stdout,*) 'TOTAL NUMBER OF ORTHONORMAL PRODUCTS:', numpw_on
     call flush_unit(stdout)

!calculate terms <\tilde{w}_j|w_i>

!     do iw=1,numpw
!        if(on_cor(iw) /= 0) then
!           omat(:,on_cor(iw))=omat(:,iw)*dsqrt(eigen(iw))
!        endif
!     enddo


     deallocate(tmpspacei,tmpspacej)

     if(.not. l_pmatrix) then
        if(ionode) then
           do iw=1,numpw_eff
              if(on_cor(iw) /=0) then
                 omat_eff(:,on_cor(iw))=omat_eff(:,iw)/dsqrt(eigen(iw))
              endif
           enddo
           allocate(omat_tmp(numpw,numpw_eff))
           omat_tmp(1:numpw,1:numpw_eff)=omat(1:numpw,numpw_offset+1:numpw_max)
           call dgemm('N','N',numpw,numpw_on, numpw_eff,1.d0,omat_tmp,numpw,omat_eff,numpw_eff,0.d0,omat,numpw)
           deallocate(omat_tmp)
        endif
        do iw=1, numpw
           call mp_barrier
           call mp_bcast(omat(1:numpw,iw),ionode_id)
        enddo
        write(stdout,*) 'OMAT 2', omat (1,1)!ATTENZIONE
     else
#ifdef __SCALAPACK
        if((myrow+1) <= nprow .and. (mycol+1) <= npcol) then
!set up omat_eff

!calculate dimesion of local omat_tmp and allocate



           numpw_on_r=ceiling(real(numpw_on)/real(nprow))
           numpw_on_c=ceiling(real(numpw_on)/real(npcol))
           numpw_on_dimr=ceiling (real(numpw_on)/real(numpw_on_r*nprow))*numpw_on_r
           numpw_on_dimc=ceiling (real(numpw_on)/real(numpw_on_c*npcol))*numpw_on_c

           write(stdout,*) 'NUMPW R',numpw, numpw_r, numpw_dimr
           write(stdout,*) 'NUMPW C',numpw, numpw_c, numpw_dimc

           write(stdout,*) 'NUMPW_EFF R',numpw_eff, numpw_eff_r, numpw_eff_dimr
           write(stdout,*) 'NUMPW_EFF C',numpw_eff, numpw_eff_c, numpw_eff_dimc

           write(stdout,*) 'NUMPW_ON R',numpw_on, numpw_on_r, numpw_on_dimr
           write(stdout,*) 'NUMPW_ON C',numpw_on, numpw_on_c, numpw_on_dimc




           allocate(omat_eff_l(numpw_eff_dimr,numpw_on_dimc))



           do iw=1,numpw_eff
              if(on_cor(iw) /=0) then
!!!!!!!!!!!!!!!!!!!!!
                 allocate(sndbuf(numpw_eff,nprow*npcol),ibuf(numpw_eff,nprow*npcol))
                 allocate(ipos(nprow*npcol))
                 ipos(:)=0
                 do jw=1,numpw_eff
                    icrow = indxg2p(jw,numpw_eff_r,0,0,nprow)
                    iccol = indxg2p(iw,numpw_eff_c,0,0,npcol)
                    iproc=icrow*npcol+iccol
                    ipos(iproc+1)=ipos(iproc+1)+1
                    ibuf(ipos(iproc+1),iproc+1)=jw
                    if(myrow==icrow .and. mycol==iccol) then
                       ilrow=indxg2l(jw,numpw_eff_r,0,0,nprow)
                       ilcol=indxg2l(iw,numpw_eff_c,0,0,npcol)
                       sndbuf(ipos(iproc+1),iproc+1)=omat_eff(ilrow,ilcol)
                    endif
                 enddo
                 do ip=1,nprow*npcol
                    if(ipos(ip)>0) then
                     !call mp_barrier
                       call mp_bcast(sndbuf(1:ipos(ip),ip),ip-1)
                       do ii=1,ipos(ip)
                          jw=ibuf(ii,ip)
                          icrow = indxg2p(jw,numpw_eff_r,0,0,nprow)
                          iccol = indxg2p(on_cor(iw),numpw_on_c,0,0,npcol)
                          if(myrow==icrow .and. mycol==iccol) then
                             ilrow=indxg2l(jw,numpw_eff_r,0,0,nprow)
                             ilcol=indxg2l(on_cor(iw),numpw_on_c,0,0,npcol)
                             omat_eff_l(ilrow,ilcol)=sndbuf(ii,ip)/dsqrt(eigen(iw))
                          endif
                       enddo
                    endif
                 enddo
                 deallocate(sndbuf,ipos,ibuf)
!!!!!!!!!!!!!!!!!!!!!!!!!
!                 do jw=1,numpw_eff
!                    icrow = indxg2p(jw,numpw_eff_r,0,0,nprow)
!                    iccol = indxg2p(iw,numpw_eff_c,0,0,npcol)
!                    iproc = icrow*npcol+iccol
!                    if(myrow==icrow .and. mycol==iccol) then
!                       ilrow=indxg2l(jw,numpw_eff_r,0,0,nprow)
!                       ilcol=indxg2l(iw,numpw_eff_c,0,0,npcol)
!                       rsca=omat_eff(ilrow,ilcol)
!                    endif
!                    call mp_bcast(rsca,iproc)
!                    icrow = indxg2p(jw,numpw_eff_r,0,0,nprow)
!                    iccol = indxg2p(on_cor(iw),numpw_on_c,0,0,npcol)
!                    if(myrow==icrow .and. mycol==iccol) then
!                       ilrow=indxg2l(jw,numpw_eff_r,0,0,nprow)
!                       ilcol=indxg2l(on_cor(iw),numpw_on_c,0,0,npcol)
!                       omat_eff_l(ilrow,ilcol)=rsca/dsqrt(eigen(iw))
!                    endif
!                 enddo
              endif
           enddo


           write(stdout,*) 'OMAT EFF L', omat_eff_l(1,1)!ATTENZIONE


!       allocate(omat_tmp(numpw,numpw_eff))

           allocate(omat_tmp(numpw_dimr, numpw_eff_dimc))
           do jw=1,numpw

              allocate(sndbuf(numpw_eff,nprow*npcol),ibuf(numpw_eff,nprow*npcol))
              allocate(ipos(nprow*npcol))
              ipos(:)=0

              do iw=1,numpw_eff
                 icrow = indxg2p(jw,numpw_r,0,0,nprow)
                 iccol = indxg2p(iw+numpw_offset,numpw_c,0,0,npcol)
                 iproc=icrow*npcol+iccol
                 ipos(iproc+1)=ipos(iproc+1)+1
                 ibuf(ipos(iproc+1),iproc+1)=iw
                 if(myrow==icrow .and. mycol==iccol) then
                    ilrow=indxg2l(jw,numpw_r,0,0,nprow)
                    ilcol=indxg2l(iw+numpw_offset,numpw_c,0,0,npcol)
                    sndbuf(ipos(iproc+1),iproc+1)=omat(ilrow,ilcol)
                 endif
              enddo
              do ip=1,nprow*npcol
                 if(ipos(ip)>0) then
                    call mp_bcast(sndbuf(1:ipos(ip),ip),ip-1)
                    do ii=1,ipos(ip)
                       iw=ibuf(ii,ip)
                       icrow = indxg2p(jw,numpw_r,0,0,nprow)
                       iccol = indxg2p(iw,numpw_eff_c,0,0,npcol)
                       if(myrow==icrow .and. mycol==iccol) then
                          ilrow=indxg2l(jw,numpw_r,0,0,nprow)
                          ilcol=indxg2l(iw,numpw_eff_c,0,0,npcol)
                          omat_tmp(ilrow,ilcol)=sndbuf(ii,ip)
                       endif
                    enddo
                 endif
              enddo
              deallocate(sndbuf,ibuf,ipos)

!                 icrow = indxg2p(jw,numpw_r,0,0,nprow)
!                 iccol = indxg2p(iw+numpw_offset,numpw_c,0,0,npcol)
!                 iproc= icrow*npcol+iccol
!                 if(myrow==icrow .and. mycol==iccol) then
!                    ilrow=indxg2l(jw,numpw_r,0,0,nprow)
!                    ilcol=indxg2l(iw+numpw_offset,numpw_c,0,0,npcol)
!                    rsca=omat(ilrow,ilcol)
!                 endif
!                 call mp_bcast(rsca,iproc)
!                 icrow = indxg2p(jw,numpw_r,0,0,nprow)
!                 iccol = indxg2p(iw,numpw_eff_c,0,0,npcol)
!                 if(myrow==icrow .and. mycol==iccol) then
!                    ilrow=indxg2l(jw,numpw_r,0,0,nprow)
!                    ilcol=indxg2l(iw,numpw_eff_c,0,0,npcol)
!                   omat_tmp(ilrow,ilcol)=rsca
!                 endif
           enddo





!pdgemm
! A = omat_tmp
             desc_a(1)=1
             desc_a(2)=icontxt
             desc_a(3)=numpw
             desc_a(4)=numpw_eff
             desc_a(5)=numpw_r
             desc_a(6)=numpw_eff_c
             desc_a(7)=0
             desc_a(8)=0
             desc_a(9)=numpw_dimr

! B = omat_eff_l

             desc_b(1)=1
             desc_b(2)=icontxt
             desc_b(3)=numpw_eff
             desc_b(4)=numpw_on
             desc_b(5)=numpw_eff_r
             desc_b(6)=numpw_on_c
             desc_b(7)=0
             desc_b(8)=0
             desc_b(9)=numpw_eff_dimr

             allocate(omat_l(numpw_dimr,numpw_on_dimc))

! C = omat_l


             desc_c(1)=1
             desc_c(2)=icontxt
             desc_c(3)=numpw
             desc_c(4)=numpw_on
             desc_c(5)=numpw_r
             desc_c(6)=numpw_on_c
             desc_c(7)=0
             desc_c(8)=0
             desc_c(9)=numpw_dimr





             call pdgemm('N','N',numpw,numpw_on,numpw_eff,1.d0,omat_tmp,1,1,desc_a,omat_eff_l,1,1,desc_b,&
                  &0.d0,omat_l,1,1,desc_c)

             !   call dgemm('N','N',numpw,numpw_on, numpw_eff,1.d0,omat_tmp,numpw,omat_eff,numpw_eff,0.d0,omat,numpw)
!collect on omat
             write(stdout,*) 'OMAT EFF 2', omat_eff(1,1)!ATTENZIONE
             write(stdout,*) 'OMAT TMP', omat_tmp(1,1)!ATTENZIONE
             write(stdout,*) 'OMAT EFF L', omat (1,1)!ATTENZIONE


             omat(:,:)=0.d0

             do jw=1,numpw
!!!!!!!!!!!!!!!!!!!!!
                allocate(sndbuf(numpw_eff,nprow*npcol),ibuf(numpw_eff,nprow*npcol))
                allocate(ipos(nprow*npcol))
                ipos(:)=0

                 do iw=1,numpw_on
                    icrow = indxg2p(jw,numpw_r,0,0,nprow)
                    iccol = indxg2p(iw,numpw_on_c,0,0,npcol)
                    iproc=icrow*npcol+iccol
                    ipos(iproc+1)=ipos(iproc+1)+1
                    ibuf(ipos(iproc+1),iproc+1)=iw
                    if(myrow==icrow .and. mycol==iccol) then
                       ilrow=indxg2l(jw,numpw_r,0,0,nprow)
                       ilcol=indxg2l(iw,numpw_on_c,0,0,npcol)
                       sndbuf(ipos(iproc+1),iproc+1)=omat_l(ilrow,ilcol)
                    endif
                 enddo
                 do ip=1,nprow*npcol
                    if(ipos(ip)>0) then
                       call mp_bcast(sndbuf(1:ipos(ip),ip),ip-1)
                       do ii=1,ipos(ip)
                          iw=ibuf(ii,ip)
                          icrow = indxg2p(jw,numpw_r,0,0,nprow)
                          iccol = indxg2p(iw,numpw_c,0,0,npcol)
                          if(myrow==icrow .and. mycol==iccol) then
                             ilrow=indxg2l(jw,numpw_r,0,0,nprow)
                             ilcol=indxg2l(iw,numpw_c,0,0,npcol)
                             omat(ilrow,ilcol)=sndbuf(ii,ip)
                          endif
                       enddo
                    endif
                 enddo
                 deallocate(sndbuf,ibuf,ipos)
!!!!!!!!!!!!!!!!!!!!!!!!!
!                do iw=1,numpw_on
!                   icrow = indxg2p(jw,numpw_r,0,0,nprow)
!                   iccol = indxg2p(iw,numpw_on_c,0,0,npcol)
!                   iproc=icrow*npcol+iccol
!                   if(myrow==icrow .and. mycol==iccol) then
!                      ilrow=indxg2l(jw,numpw_r,0,0,nprow)
!                      ilcol=indxg2l(iw,numpw_on_c,0,0,npcol)
!                      rsca=omat_l(ilrow,ilcol)
!                   endif
!                   call mp_barrier!ATTENZIONE
!                   call mp_bcast(rsca,iproc)
!                    icrow = indxg2p(jw,numpw_r,0,0,nprow)
!                    iccol = indxg2p(iw,numpw_c,0,0,npcol)
!                    if(myrow==icrow .and. mycol==iccol) then
!                       ilrow=indxg2l(jw,numpw_r,0,0,nprow)
!                       ilcol=indxg2l(iw,numpw_c,0,0,npcol)
!                       omat(ilrow,ilcol)=rsca
!                    endif
!                enddo
              enddo

             write(stdout,*) 'OMAT 2', omat (1,1)!ATTENZIONE
             call flush_unit(stdout)

           deallocate(omat_eff_l,omat_l)

           deallocate(omat_tmp)
        else

        endif

#endif
     endif

     numpw_original=numpw
     numpw=numpw_on
!ATTENTION the dimension of the omat matrix remains UNCHANGED
     if(l_pmatrix) then
        deallocate (omat_eff)
     else
        if(ionode) deallocate(omat_eff)
        deallocate(omat_eff_bcast)
     endif
     deallocate(eigen, on_cor)
     return
   end subroutine orthonormalize_producs_cutoff
  !
  ! -------------------------------------------------
  !
   subroutine do_polarization_analysis( wp, num_wp_dimr, num_wp_r,mcut,numpw,numpw_original, &
                                n_set, ecutoff, lcutoff,omat,newdimr, newdimc,nc_max,n_r,n_c)
!this subroutine
!1)construct the polarization matrix at zero time
!2)find the  eigenvalues above cutoff
!3)redefines the products of wannier functions
!it also updats the .gzero file
!it also deallocate the descriptors wp

     USE kinds,      ONLY : DP
     USE wannier_gw, ONLY : wannier_P, u_trans,num_nbndv, free_memory,nbnd_normal,max_ngm, &
          &l_pmatrix,npcol,nprow,myrow,mycol,icontxt,l_assume_ortho,l_coulomb_analysis
 !u_trans is defined as Psi_i=U_{i,j}w_j
     USE io_global, ONLY : stdout, ionode, ionode_id
     USE io_files,  ONLY : find_free_unit, prefix, diropn
     USE gvect,     ONLY : ngm, gg,gstart
     USE cell_base, ONLY: tpiba2
     USE mp,         ONLY : mp_bcast, mp_sum, mp_barrier
     USE mp_global,  ONLY : mpime, nproc

     implicit none

     INTEGER, INTENT(inout) :: numpw!dimension of matrices in  MODIFIED IN OUTPUT
     INTEGER, INTENT(in) :: numpw_original!dimension of original matrices
     INTEGER :: n_set !length of buffer for products of wannier
     REAL(kind=DP), INTENT(in) :: ecutoff!cutoff in Rydberg for g sum
     LOGICAL, INTENT(in) :: lcutoff !if true uses cutoff on G defined by ecutoff POTENTIALLY DANGEROUS
     TYPE(wannier_P), INTENT(inout) :: wp(num_wp_dimr)!descriptor of products of wanniers
     INTEGER, INTENT(in) :: num_wp_dimr!dimension of wp
     INTEGER, INTENT(in) :: num_wp_r!periodicity of wp for parallel matrix execution
     REAL(kind=DP), INTENT(in)      :: mcut!cutoff for eigenvalues of the polarization matrix
     INTEGER, INTENT(in) :: newdimr,newdimc!dimensions of omat
     REAL(kind=DP), INTENT(inout) :: omat(newdimr,newdimc)!overlap matrix <\tilde{w}^P_i|{w}^P'_j>
                                                    !to be updated
     INTEGER, INTENT(in) :: nc_max!number of conduction states to be used for building the polarization matrix
     INTEGER, INTENT(in) :: n_r, n_c!in parallel case periodicity on rows and columns


!internal variables

     REAL(kind=DP), ALLOCATABLE :: contv_ci(:,:), pola(:,:)
     INTEGER :: nbndc
     INTEGER :: iv,ic,iw,jw, ii, jj,kk,ll,ip
     REAL(kind=DP), ALLOCATABLE :: eigen(:),work(:)
     INTEGER :: lwork,info

     INTEGER :: i,j,k, ig
     REAL(kind=DP) :: fact
     INTEGER :: iunu, iungprod,iungprod_red
     COMPLEX(kind=DP), ALLOCATABLE :: tmpspacei(:,:)
     COMPLEX(kind=DP), ALLOCATABLE :: tmpspacej(:,:)
     INTEGER :: iiw,jjw
     INTEGER, ALLOCATABLE :: on_cor(:)
     INTEGER :: numpw_red
     REAL(kind=DP), ALLOCATABLE :: tmp_str(:)
     LOGICAL :: exst
     INTEGER :: ngm_max
     INTEGER :: iun_g
     REAL(kind=DP), ALLOCATABLE :: gzero(:), omat_tmp(:,:), pola_tmp(:,:)
     INTEGER :: numpw_r,numpw_c,numpw_dimr,numpw_dimc,numpw_original_r,numpw_original_c,numpw_original_dimr,numpw_original_dimc
     INTEGER :: icrow,iccol,ilrow,ilcol,iproc
     REAL(kind=DP) :: rsca
#ifdef __SCALAPACK
     INTEGER, EXTERNAL :: indxg2p,indxg2l
#endif
     INTEGER :: desc_a(9),desc_b(9),desc_c(9)

     REAL(kind=DP), ALLOCATABLE :: sndbuf(:,:)
     INTEGER, ALLOCATABLE :: ipos(:),ibuf(:,:)


     write(stdout,*) 'NCMAX',nc_max!ATTENZIONE
     write(stdout,*) 'DO POLARIZATION ANALYSIS:',numpw,numpw_original
     write(stdout,*) 'DO POLARIZATION ANALYSIS:',newdimr, newdimc,nc_max,n_r,n_c
     call flush_unit(stdout)

!allocate and set to zero
     nbndc=nbnd_normal-num_nbndv
     allocate(pola(numpw,numpw))
     allocate(contv_ci(nbndc,numpw))
     allocate(gzero(numpw))
     gzero(:)=0.d0
     pola(:,:)=0.d0

!loop on v
     if(.not.l_assume_ortho) then
        do iv=1,num_nbndv
!        if(mod(iv,nproc)==mpime) then
           contv_ci(:,:)=0.d0
!calculate contraction
           do iw=1,numpw
              if(.not.l_pmatrix) then
                 if(ionode) then
                    do jw=1,wp(iw)%numij
                       ii=wp(iw)%ij(1,jw)
                       jj=wp(iw)%ij(2,jw)
                       if(ii>num_nbndv) then
                          write(stdout,*) 'ERRORE ON ORDER WP'
                          call flush_unit(stdout)
                          stop
                       endif
                       if(jj>num_nbndv) then!IMPORTANT the definiton of product order
                          do ll=1,min(nbndc,nc_max)
                             contv_ci(ll,iw)=contv_ci(ll,iw)+ &
                                  &dble(u_trans(iv,ii))*dble(u_trans(ll+num_nbndv,jj))*wp(iw)%o(jw)
                          enddo
                       endif
                    enddo
                 endif
                 call mp_bcast(contv_ci(:,iw),ionode_id)
              else
#ifdef __SCALAPACK
                 icrow=indxg2p(iw,num_wp_r,0,0,nproc)
                 if(icrow==mpime) then
                    ilrow=indxg2l(iw,num_wp_r,0,0,nproc)
                    do jw=1,wp(ilrow)%numij
                       ii=wp(ilrow)%ij(1,jw)
                       jj=wp(ilrow)%ij(2,jw)
                       if(ii>num_nbndv) then
                          write(stdout,*) 'ERRORE ON ORDER WP'
                          call flush_unit(stdout)
                          stop
                       endif
                       if(jj>num_nbndv) then!IMPORTANT the definiton of product order
                          do ll=1,min(nbndc,nc_max)
                             !   if(jj<=num_nbndv+min(nbndc,nc_max)) then
                             contv_ci(ll,iw)=contv_ci(ll,iw)+ &
                                  &dble(u_trans(iv,ii))*dble(u_trans(ll+num_nbndv,jj))*wp(ilrow)%o(jw)
                             !  endif
                          enddo
                       endif
                    enddo
                 endif
#endif
              endif
           enddo
           if(l_pmatrix) call mp_sum(contv_ci(:,:))
!loop on i,j,c
!calculate polarization
           do iw=1,numpw
              do jw=1,numpw
                 do ic=1,min(nbndc,nc_max)
                    pola(jw,iw)=pola(jw,iw)+contv_ci(ic,iw)*contv_ci(ic,jw)
                 enddo
              enddo
           enddo
!        endif
        enddo
     endif
!     do iw=1,numpw
!        call mp_sum(pola(:,iw))
!     enddo
     deallocate(contv_ci)

!calculate eigen values eigen vectors of polarization matrix
!apply cut off
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     write(stdout,*) 'MATRIX PREPARED'
     call flush_unit(stdout)


     if(lcutoff) then
        ngm_max=0
        do ig=1,ngm
           if(gg(ig)*tpiba2 >= ecutoff) exit
           ngm_max=ngm_max+1
        enddo
     else
        ngm_max=ngm
     endif

     write(stdout,*) 'NGM MAX:', ngm_max, ngm






     allocate(eigen(numpw))
     if(.not. l_assume_ortho) then
        if(ionode) then
           allocate(work(1))
           call  DSYEV( 'V', 'U', numpw, pola, numpw, eigen, work, -1, info )
           lwork=work(1)
           deallocate(work)
           allocate(work(lwork))
           call  DSYEV( 'V', 'U', numpw, pola, numpw, eigen, work, lwork, info )
           deallocate(work)
           if(info/=0) then
              write(stdout,*) 'ROUTINE find_transform, INFO:', info
              stop
           endif
        endif
        do iw=1,numpw
           call mp_barrier
           call mp_bcast(pola(:,iw), ionode_id)
        enddo
        call mp_bcast(eigen(:), ionode_id)
     else
        pola(:,:)=0.d0
        do i=1,numpw
           pola(i,i)=1.d0
        enddo
        do i=1,numpw
           eigen(i)=100.d0
        enddo
     endif

     do iw=1,numpw
        write(stdout,*) 'EIGEN:',iw, eigen(iw)
     enddo



     iungprod = find_free_unit()
     CALL diropn( iungprod, 'wiwjwfc_on', max_ngm*2, exst )
     iungprod_red = find_free_unit()
     if(.not.l_coulomb_analysis) then
        CALL diropn( iungprod_red, 'wiwjwfc_red', max_ngm*2, exst )
     else
        CALL diropn( iungprod_red, 'wiwjwfc_red_1', max_ngm*2, exst )
     endif

     if(ionode) then
        iun_g =  find_free_unit()
        open( unit= iun_g, file=trim(prefix)//'.gzero', status='unknown',form='unformatted')
     endif




!allocate file
     allocate(tmpspacei(max_ngm,n_set),tmpspacej(max_ngm,n_set))
     allocate(on_cor(numpw))
     on_cor(:)=0

     numpw_red=0
!external loop on orthonormal products
     do iiw=1,ceiling(real(numpw)/real(n_set))
!if required read them from disk or set to zero
        write(stdout,*) 'IIW', iiw
        call flush_unit(stdout)
        tmpspacei(:,:)=0.d0
        do jjw=1,ceiling(real(numpw)/real(n_set))
            write(stdout,*) 'JJW', jjw
            call flush_unit(stdout)
            do jw=(jjw-1)*n_set+1,min(jjw*n_set,numpw)
               CALL davcio(tmpspacej(:,jw-(jjw-1)*n_set), max_ngm*2,iungprod,jw,-1)
            enddo
            do iw=(iiw-1)*n_set+1,min(iiw*n_set,numpw)
               if(eigen(iw) >= mcut) then
                  if(on_cor(iw) == 0) then
                     numpw_red=numpw_red+1
                     on_cor(iw)=numpw_red
                  endif
                  do jw=(jjw-1)*n_set+1,min(jjw*n_set,numpw)
                     tmpspacei(1:ngm_max,iw-(iiw-1)*n_set)=tmpspacei(1:ngm_max,iw-(iiw-1)*n_set)+pola(jw,iw)*&
                          &tmpspacej(1:ngm_max,jw-(jjw-1)*n_set)
                  enddo
               endif
            enddo
         enddo
         do iw=(iiw-1)*n_set+1,min(iiw*n_set,numpw)
            if(on_cor(iw) /= 0) then
               CALL davcio(tmpspacei(:,iw-(iiw-1)*n_set), max_ngm*2,iungprod_red,on_cor(iw),1)
               if(gstart==2) gzero(on_cor(iw))=dble(tmpspacei(1,iw-(iiw-1)*n_set))
            endif
         enddo
      enddo



     close(iungprod)
     close(iungprod_red)


!update gzero file
     call mp_sum(gzero(:))
     if(ionode) then
        write(iun_g) numpw_red
        do iw=1,numpw_red
           write(iun_g) gzero(iw)
           if(ionode) write(stdout,*) 'GZERO:',iw,gzero(iw),numpw_red
           call flush_unit(stdout)
        enddo
        close(iun_g)
     endif
!for compatibility write orthonorm file


     if(ionode) then
        allocate(tmp_str(numpw_red))
        iunu = find_free_unit()

        open(unit=iunu,file=trim(prefix)//'.orthonorm',status='unknown',form='unformatted')
        write(iunu) numpw_red

        do  i=1,numpw_red
           tmp_str(:)=0.d0
           tmp_str(i)=1.d0
           write(iunu) tmp_str
        enddo

        close(iunu)
        deallocate(tmp_str)
     endif

     write(stdout,*) 'TOTAL NUMBER OF ORTHONORMAL PRODUCTS:', numpw_red
     call flush_unit(stdout)

!update overlap matrix omat
     deallocate(tmpspacei,tmpspacej)

!deallocate descriptors of wannier products

     do iw=1,numpw
        if(.not.l_pmatrix) then
           if(ionode) call free_memory(wp(iw))
        else
#ifdef __SCALAPACK
           icrow=indxg2p(iw,num_wp_r,0,0,nproc)
           if(icrow==mpime) then
              ilrow=indxg2l(iw,num_wp_r,0,0,nproc)
              call free_memory(wp(ilrow))
           endif
#endif
        endif
     enddo


     if(.not.l_pmatrix) then
        if(ionode)allocate(omat_tmp(numpw_original,numpw))
     else
        numpw_r=n_r!ceiling(real(numpw)/real(max(nprow,npcol)))
        numpw_c=n_c!ceiling(real(numpw)/real(max(nprow,npcol)))
        numpw_dimr=ceiling (real(numpw)/real(numpw_r*nprow))*numpw_r
        numpw_dimc=ceiling (real(numpw)/real(numpw_c*npcol))*numpw_c
        numpw_original_r=n_r!ceiling(real(numpw_original)/real(nprow))
        numpw_original_c=n_c!ceiling(real(numpw_original)/real(npcol))
        numpw_original_dimr=newdimr!ceiling (real(numpw_original)/real(numpw_original_r*nprow))*numpw_original_r
        numpw_original_dimc=newdimc!ceiling (real(numpw_original)/real(numpw_original_c*npcol))*numpw_original_c
       allocate(omat_tmp(numpw_original_dimr,numpw_dimc))
     endif
     do iw=1,numpw
        if(on_cor(iw)/=0) then
           pola(:,on_cor(iw))=pola(:,iw)
        endif
     enddo

     if(.not.l_pmatrix) then
        if(ionode) omat_tmp(1:numpw_original, 1:numpw)=omat(1:numpw_original,1:numpw)
     else
#ifdef __SCALAPACK
        do jw=1,numpw_original

           allocate(sndbuf(numpw,nprow*npcol),ibuf(numpw,nprow*npcol))
           allocate(ipos(nprow*npcol))
           ipos(:)=0
           do iw=1,numpw
              icrow = indxg2p(jw,n_r,0,0,nprow)
              iccol = indxg2p(iw,n_c,0,0,npcol)
              iproc=icrow*npcol+iccol
              ipos(iproc+1)=ipos(iproc+1)+1
              ibuf(ipos(iproc+1),iproc+1)=iw
              if(myrow==icrow .and. mycol==iccol) then
                 ilrow=indxg2l(jw,n_r,0,0,nprow)
                 ilcol=indxg2l(iw,n_c,0,0,npcol)
                 sndbuf(ipos(iproc+1),iproc+1)=omat(ilrow,ilcol)
              endif
           enddo
           do ip=1,nprow*npcol
              if(ipos(ip)>0) then
                 call mp_bcast(sndbuf(1:ipos(ip),ip),ip-1)
                 do ii=1,ipos(ip)
                    iw=ibuf(ii,ip)
                    icrow = indxg2p(jw,numpw_original_r,0,0,nprow)
                    iccol = indxg2p(iw,numpw_c,0,0,npcol)
                    if(myrow==icrow .and. mycol==iccol) then
                       ilrow=indxg2l(jw,numpw_original_r,0,0,nprow)
                       ilcol=indxg2l(iw,numpw_c,0,0,npcol)
                       omat_tmp(ilrow,ilcol)=sndbuf(ii,ip)
                    endif
                 enddo
              endif
           enddo
           deallocate(sndbuf,ibuf,ipos)
!
!           do iw=1,numpw
!              icrow = indxg2p(jw,n_r,0,0,nprow)
!              iccol = indxg2p(iw,n_c,0,0,npcol)
!              iproc=icrow*npcol+iccol
!              if(myrow==icrow .and. mycol==iccol) then
!                 ilrow=indxg2l(jw,n_r,0,0,nprow)
!                 ilcol=indxg2l(iw,n_c,0,0,npcol)
!                 rsca=omat(ilrow,ilcol)
!              endif
!              call mp_bcast(rsca,iproc)
!              icrow = indxg2p(jw,numpw_original_r,0,0,nprow)
!              iccol = indxg2p(iw,numpw_c,0,0,npcol)
!              if(myrow==icrow .and. mycol==iccol) then
!                 ilrow=indxg2l(jw,numpw_original_r,0,0,nprow)
!                 ilcol=indxg2l(iw,numpw_c,0,0,npcol)
!                 omat_tmp(ilrow,ilcol)=rsca
!              endif
!           enddo
        enddo
#endif
     endif

     if(.not.l_pmatrix) then
        if(ionode)call dgemm('N','N',numpw_original,numpw_red, numpw,1.d0,omat_tmp,&
          numpw_original,pola,numpw,0.d0,omat,numpw_original)
        do iw=1,numpw_red
           call mp_barrier
           call mp_bcast(omat(1:numpw_original,iw),ionode_id)
        enddo
     else
#ifdef __SCALAPACK
!fatto fin qui
        allocate(pola_tmp(numpw_dimr,newdimc))
        do jw=1,numpw
           do iw=1,numpw
              icrow = indxg2p(jw,numpw_r,0,0,nprow)
              iccol = indxg2p(iw,n_c,0,0,npcol)
              if(icrow == myrow .and. iccol == mycol) then
                 ilrow = indxg2l(jw,numpw_r,0,0,nprow)
                 ilcol = indxg2l(iw,n_c,0,0,npcol)
                 pola_tmp(ilrow,ilcol)=pola(jw,iw)
              endif
           enddo
        enddo

!pdgemm
! A = omat_tmp
             desc_a(1)=1
             desc_a(2)=icontxt
             desc_a(3)=numpw_original
             desc_a(4)=numpw
             desc_a(5)=numpw_original_r
             desc_a(6)=numpw_c
             desc_a(7)=0
             desc_a(8)=0
             desc_a(9)=numpw_original_dimr

! B = pola_tmp
             desc_b(1)=1
             desc_b(2)=icontxt
             desc_b(3)=numpw
             desc_b(4)=numpw
             desc_b(5)=numpw_r
             desc_b(6)=n_c
             desc_b(7)=0
             desc_b(8)=0
             desc_b(9)=numpw_dimr


! C = omat
             desc_c(1)=1
             desc_c(2)=icontxt
             desc_c(3)=numpw_original
             desc_c(4)=numpw_red
             desc_c(5)=n_r
             desc_c(6)=n_c
             desc_c(7)=0
             desc_c(8)=0
             desc_c(9)=newdimr

             call pdgemm('N','N',numpw_original,numpw_red, numpw,1.d0,omat_tmp,1,1,desc_a,pola_tmp,1,1,desc_b,0.d0,omat,&
                  &1,1,desc_c)


        deallocate(pola_tmp)
#endif
     endif
     if(ionode) deallocate(omat_tmp)

!update dimension of reduced orthonormal products of wanniers
     numpw=numpw_red








     deallocate(eigen, on_cor)
     deallocate(gzero)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     deallocate(pola)
     write(stdout,*) 'EXITING do_polarization_analysis'
     call flush_unit(stdout)
  return
end subroutine do_polarization_analysis
 !
 ! ------------------------------------------------
 !
subroutine do_polarization_analysis_study( wp,mcut,numpw,numpw_original, n_set, ecutoff, lcutoff,omat,nc_max)
!this subroutine
!1)construct the polarization matrix at zero time
!2)find the  eigenvalues above cutoff
!3)redefines the products of wannier functions
!it also updats the .gzero file
!it also deallocate the descriptors wp

     USE kinds,      ONLY : DP
     USE wannier_gw, ONLY : wannier_P, u_trans,num_nbndv, free_memory,max_ngm
     USE wvfct,      ONLY : nbnd, et
 !u_trans is defined as Psi_i=U_{i,j}w_j
     USE io_global, ONLY : stdout, ionode, ionode_id
     USE io_files,  ONLY : find_free_unit, prefix
     USE gvect,     ONLY : ngm, gg,gstart
     USE cell_base, ONLY: tpiba2
     USE mp,         ONLY : mp_bcast, mp_sum, mp_barrier
     USE mp_global,  ONLY : mpime, nproc

     implicit none

     INTEGER, INTENT(inout) :: numpw!dimension of matrices in  MODIFIED IN OUTPUT
     INTEGER, INTENT(in) :: numpw_original!dimension of original matrices
     INTEGER :: n_set !length of buffer for products of wannier
     REAL(kind=DP), INTENT(in) :: ecutoff!cutoff in Rydberg for g sum
     LOGICAL, INTENT(in) :: lcutoff !if true uses cutoff on G defined by ecutoff POTENTIALLY DANGEROUS
     TYPE(wannier_P), INTENT(inout) :: wp(numpw)!descriptor of products of wanniers
     REAL(kind=DP), INTENT(in)      :: mcut!cutoff for eigenvalues of the polarization matrix
     REAL(kind=DP), INTENT(inout) :: omat(numpw_original, numpw_original)!overlap matrix <\tilde{w}^P_i|{w}^P'_j>
                                                    !to be updated
     INTEGER, INTENT(in) :: nc_max!number of conduction states to be used for building the polarization matrix


!internal variables

     REAL(kind=DP), ALLOCATABLE :: contv_ci(:,:), pola(:,:,:)
     INTEGER :: nbndc
     INTEGER :: iv,ic,iw,jw, ii, jj,kk,ll, it, ie
     REAL(kind=DP), ALLOCATABLE :: eigen(:),work(:)
     INTEGER :: lwork,info

     INTEGER :: i,j,k, ig
     REAL(kind=DP) :: fact
     INTEGER :: iunu, iungprod,iungprod_red
     COMPLEX(kind=DP), ALLOCATABLE :: tmpspacei(:,:)
     COMPLEX(kind=DP), ALLOCATABLE :: tmpspacej(:,:)
     INTEGER :: iiw,jjw
     INTEGER :: numpw_red
     REAL(kind=DP), ALLOCATABLE :: tmp_str(:)
     LOGICAL :: exst
     INTEGER :: ngm_max
     INTEGER :: iun_g
     REAL(kind=DP), ALLOCATABLE :: o_proj(:,:,:), o_proj_mul(:,:), o_proj_tmp(:,:)

     INTEGER :: iun_proj
     INTEGER :: n_times = 10
     INTEGER :: n_first = 250
     REAL(kind=DP) :: times(10), trace

     times(1)=0.d0
     times(2)=1.d0
     times(3)=2.d0
     times(4)=3.d0
     times(5)=4.d0
     times(6)=5.d0
     times(7)=6.d0
     times(8)=7.d0
     times(9)=8.d0
     times(10)=9.d0

     write(stdout,*) 'NCMAX',nc_max!ATTENZIONE

     if(ionode) then
        iun_proj =  find_free_unit()
        open(unit=iun_proj,file='projections.dat',status='unknown')
     endif

!allocate and set to zero
     nbndc=nbnd-num_nbndv
     allocate(pola(numpw,numpw,n_times))
     allocate(contv_ci(nbndc,numpw))
     allocate(o_proj(numpw,numpw,n_times))

     do it=1,n_times
        pola(:,:,it)=0.d0

!loop on v
        do iv=1,num_nbndv
           if(mod(iv,nproc)==mpime) then
              contv_ci(:,:)=0.d0
!calculate contraction
              do iw=1,numpw
                 if(ionode) then
                    do jw=1,wp(iw)%numij
                       ii=wp(iw)%ij(1,jw)
                       jj=wp(iw)%ij(2,jw)
                       if(ii>num_nbndv) then
                          write(stdout,*) 'ERRORE ON ORDER WP'
                          call flush_unit(stdout)
                          stop
                       endif
                       if(jj>num_nbndv) then!IMPORTANT the definiton of product order
                          do ll=1,min(nbndc,nc_max)
                             fact=dsqrt(exp(-(et(ll+num_nbndv,1)-et(iv,1))*abs(times(it))))
                     !   if(jj<=num_nbndv+min(nbndc,nc_max)) then
                             contv_ci(ll,iw)=contv_ci(ll,iw)+ &
                                  &dble(u_trans(iv,ii))*dble(u_trans(ll+num_nbndv,jj))*wp(iw)%o(jw)*fact
                     !  endif
                          enddo
                       endif
                    enddo
                 endif
                 call mp_bcast(contv_ci(:,iw),ionode_id)
              enddo
!loop on i,j,c
!calculate polarization
              do iw=1,numpw
                 do jw=1,numpw
                    do ic=1,min(nbndc,nc_max)
                       pola(jw,iw,it)=pola(jw,iw,it)+contv_ci(ic,iw)*contv_ci(ic,jw)
                    enddo
                 enddo
              enddo
           endif
        enddo
        do iw=1,numpw
           call mp_sum(pola(:,iw,it))
        enddo


!calculate eigen values eigen vectors of polarization matrix
!apply cut off
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!






        allocate(eigen(numpw))
        if(ionode) then
           allocate(work(1))
           call  DSYEV( 'V', 'U', numpw, pola(:,:,it), numpw, eigen, work, -1, info )
           lwork=work(1)
           deallocate(work)
           allocate(work(lwork))
           call  DSYEV( 'V', 'U', numpw, pola(:,:,it), numpw, eigen, work, lwork, info )
           deallocate(work)
           if(info/=0) then
              write(stdout,*) 'ROUTINE find_transform, INFO:', info
              stop
           endif
        endif
        do iw=1,numpw
           call mp_bcast(pola(:,iw,it), ionode_id)
        enddo
        call mp_bcast(eigen(:), ionode_id)

        do iw=1,numpw
           write(stdout,*) 'EIGEN:',iw, eigen(iw)
           call flush_unit(stdout)
        enddo

        o_proj(:,:,it)=0.d0

        deallocate(eigen)
        if(ionode) then
           write(iun_proj,*) it
           do iw=1,numpw
              do jw=1,numpw
                 write(iun_proj,*) pola(jw,iw,it)
              enddo
           enddo
        endif
     enddo

     if(ionode) close(iun_proj)
     do ie=1,n_first
        do it=1,n_times
           o_proj(:,:,it)=0.d0
           do ii=1,ie
              do iw=1,numpw
                 do jw=1,numpw
                    o_proj(iw,jw,it)=o_proj(iw,jw,it)+pola(iw,numpw-ii+1,it)*pola(jw,numpw-ii+1,it)
                 enddo
              enddo
           enddo


        enddo
        allocate(o_proj_mul(numpw,numpw))
        allocate(o_proj_tmp(numpw,numpw))
        call dgemm('N','N',numpw,numpw,numpw,1.d0,o_proj(:,:,1),numpw,o_proj(:,:,2),numpw,0.d0,o_proj_mul,numpw)
!        call dgemm('N','N',numpw,numpw,numpw,1.d0,o_proj_tmp,numpw,o_proj(:,:,3),numpw,0.d0,o_proj_mul,numpw)

        trace=0.d0
        do iw=1,numpw
           trace=trace+o_proj_mul(iw,iw)
        enddo

        write(stdout,*) 'TRACE:',ie, trace
           !deallocate descriptors of wannier products

        deallocate(o_proj_mul,o_proj_tmp)
     enddo

     deallocate(contv_ci)
     deallocate(o_proj)
     deallocate(pola)

  return
end subroutine do_polarization_analysis_study
 !
 ! ------------------------------------------------
 !
  subroutine coulomb_analysis(numpw,numpw_original,numpw_vvc,vmat, numpwdimr,numpwdimc,cutoff, n_set, ecutoff, lcutoff,&
                                  &omat, newdimr,newdimc, n_r, n_c)

!this subroutine calculate the eigenvalues and eigenvectors of the coulomb
!potential vmat and then only the states corresponding to larger eigenvalues are taken

     USE kinds, ONLY : DP
     USE io_global, ONLY : stdout, ionode, ionode_id
     USE io_files, ONLY : find_free_unit, prefix, diropn
     USE gvect,  ONLY : ngm, gg
     USE cell_base, ONLY: tpiba2
     USE mp,        ONLY : mp_bcast, mp_sum,mp_barrier
     USE wannier_gw, ONLY : max_ngm,p_mpime,p_nproc, npcol, nprow,icontxt,myrow,mycol,l_pmatrix

     implicit none

     INTEGER, INTENT(inout) :: numpw!dimension of matrices in  MODIFIED IN OUTPUT
     INTEGER, INTENT(out) :: numpw_original!dimension of original matrices
     REAL(kind=DP) :: vmat(numpwdimr,numpwdimc)!overlap matrix ON OUTPUT contains the terms <\tilde{w}_j|w_i>
     INTEGER,INTENT(in) :: numpwdimr,numpwdimc!dimensions of omat
     REAL(kind=DP) :: cutoff!cutoff for eigenvalues
     INTEGER :: n_set !length of buffer for products of wannier
     REAL(kind=DP), INTENT(in) :: ecutoff!cutoff in Rydberg for g sum
     LOGICAL, INTENT(in) :: lcutoff !if true uses cutoff on G defined by ecutoff POTENTIALLY DANGEROUS
     REAL(kind=DP), INTENT(inout) :: omat(newdimr,newdimc)!overlap matrix to be updated
     INTEGER, INTENT(in) :: newdimr,newdimc!dimensions of omat
     INTEGER, INTENT(in) :: n_r, n_c!periodicity of omat
     INTEGER :: numpw_vvc!total number of non-orthonomal wannier products


     REAL(kind=DP), ALLOCATABLE :: eigen(:),work(:)
     INTEGER :: lwork,info,liwork
     INTEGER, ALLOCATABLE :: iwork(:)

     INTEGER :: i,j,k, ig
     REAL(kind=DP) :: fact
     INTEGER :: iunu, iungprod,iungprod_on
     COMPLEX(kind=DP), ALLOCATABLE :: tmpspacei(:,:)
     COMPLEX(kind=DP), ALLOCATABLE :: tmpspacej(:,:)
     INTEGER :: iw,jw,iiw,jjw
     INTEGER, ALLOCATABLE :: on_cor(:)
     INTEGER :: numpw_on
     REAL(kind=DP), ALLOCATABLE :: tmp_str(:)
     LOGICAL :: exst
     INTEGER :: ngm_max




     INTEGER :: numpw_r,numpw_c,numpw_dimr,numpw_dimc
     INTEGER :: numpw_vvc_r,numpw_vvc_c,numpw_vvc_dimr,numpw_vvc_dimc
     REAL(kind=DP), ALLOCATABLE :: vmat_states(:,:), omat_tmp(:,:)
     INTEGER :: desc_a(9),desc_b(9),desc_c(9)
     INTEGER, ALLOCATABLE :: ifail(:),iclustr(:)
     REAL(kind=DP), ALLOCATABLE :: gap(:)
     INTEGER :: m,nz,icrow,iccol,iproc,ilrow,ilcol
#ifdef __SCALAPACK
     INTEGER, EXTERNAL :: indxg2p,indxg2l
#endif
     REAL(kind=DP) :: rsca

     write(stdout,*) 'ROUTINE coulomb_analysis'



     !determine ngm_max
     if(lcutoff) then
        ngm_max=0
        do ig=1,ngm
           if(gg(ig)*tpiba2 >= ecutoff) exit
           ngm_max=ngm_max+1
        enddo
     else
        ngm_max=ngm
     endif

     write(stdout,*) 'NGM MAX:', ngm_max, ngm
     call flush_unit(stdout)


     if(.not.l_pmatrix) then





        allocate(eigen(numpw))
        if(ionode) then
           allocate(work(1))
           call  DSYEV( 'V', 'U', numpw, vmat, numpwdimr, eigen, work, -1, info )
           lwork=work(1)
           deallocate(work)
           allocate(work(lwork))
           call  DSYEV( 'V', 'U', numpw, vmat, numpwdimr, eigen, work, lwork, info )
           deallocate(work)
           if(info/=0) then
              write(stdout,*) 'ROUTINE find_transform, INFO:', info
              stop
           endif
        endif
        do iw=1,numpw
           call mp_bcast(vmat(:,iw), ionode_id)
        enddo
        call mp_bcast(eigen(:), ionode_id)
     else
#ifdef __SCALAPACK



        numpw_r=ceiling(real(numpw)/real(nprow))
        numpw_c=ceiling(real(numpw)/real(npcol))
        numpw_dimr=ceiling (real(numpw)/real(numpw_r*nprow))*numpw_r
        numpw_dimc=ceiling (real(numpw)/real(numpw_c*npcol))*numpw_c

        if((myrow+1) <= nprow .and. (mycol+1) <= npcol) then

           allocate(eigen(numpw))
           allocate(gap(numpw),ifail(numpw),iclustr(numpw))
           allocate(vmat_states(numpw_dimr,numpw_dimc))





! A = vmat
             desc_a(1)=1
             desc_a(2)=icontxt
             desc_a(3)=numpw
             desc_a(4)=numpw
             desc_a(5)=numpw_r
             desc_a(6)=numpw_c
             desc_a(7)=0
             desc_a(8)=0
             desc_a(9)=numpw_dimr
!B= vmat_states
             desc_b(1)=1
             desc_b(2)=icontxt
             desc_b(3)=numpw
             desc_b(4)=numpw
             desc_b(5)=numpw_r
             desc_b(6)=numpw_c
             desc_b(7)=0
             desc_b(8)=0
             desc_b(9)=numpw_dimr

             allocate(work(1),iwork(1))




             call pdsyevx('V','A','U',numpw,vmat,1,1,desc_a,0.d0,0.d0,0,0,0.d0,&
                  & m,nz,eigen,-1.d0,vmat_states,1,1,desc_b,work,-1,iwork,-1,ifail,iclustr,gap,info)




             lwork=work(1)
             liwork=iwork(1)
             deallocate(work,iwork)
             write(stdout,*) 'LWORK LIWORK:',lwork,liwork
             call flush_unit(stdout)
             allocate(work(lwork),iwork(liwork))
             call pdsyevx('V','A','U',numpw,vmat,1,1,desc_a,0.d0,0.d0,0,0,0.d0,&
       & m,nz,eigen,0.d0,vmat_states,1,1,desc_b,work,lwork,iwork,liwork,ifail,iclustr,gap,info)
             if(info/=0) then
                write(stdout,*) 'ROUTINE PDSYEV INFO:', info
                stop
             endif


             deallocate(work,iwork)
             deallocate(gap,iclustr,ifail)


          endif
#endif
       endif



       call mp_bcast(eigen(:),ionode_id)


       do iw=1,numpw
          write(stdout,*) 'EIGEN:',iw, eigen(iw)
       enddo
       call flush_unit(stdout)


       iungprod = find_free_unit()
       CALL diropn( iungprod, 'wiwjwfc_red_1', max_ngm*2, exst )
       iungprod_on = find_free_unit()
       CALL diropn( iungprod_on, 'wiwjwfc_red', max_ngm*2, exst )



!allocate file
       allocate(tmpspacei(max_ngm,n_set),tmpspacej(max_ngm,n_set))
       allocate(on_cor(numpw))
       on_cor(:)=0

       numpw_on=0
!external loop on orthonormal products
       do iiw=1,ceiling(real(numpw)/real(n_set))
!if required read them from disk or set to zero
          write(stdout,*) 'IIW', iiw
          call flush_unit(stdout)
          tmpspacei(:,:)=0.d0
          do jjw=1,ceiling(real(numpw)/real(n_set))
             write(stdout,*) 'JJW', jjw
             call flush_unit(stdout)
             do jw=(jjw-1)*n_set+1,min(jjw*n_set,numpw)
                CALL davcio(tmpspacej(:,jw-(jjw-1)*n_set), max_ngm*2,iungprod,jw,-1)
             enddo
             do iw=(iiw-1)*n_set+1,min(iiw*n_set,numpw)
                if(eigen(iw) >= cutoff) then
                   if(on_cor(iw) == 0) then
                      numpw_on=numpw_on+1
                      on_cor(iw)=numpw_on
                   endif
                   do jw=(jjw-1)*n_set+1,min(jjw*n_set,numpw)
                      if(.not.l_pmatrix) then
                         tmpspacei(1:ngm_max,iw-(iiw-1)*n_set)=tmpspacei(1:ngm_max,iw-(iiw-1)*n_set)+vmat(jw,iw)*&
                              &tmpspacej(1:ngm_max,jw-(jjw-1)*n_set)
                      else
#ifdef __SCALAPACK
!determine if jw and iw stays in the processor
!USARE LAPACK ROUTINES
                         icrow = indxg2p(jw,numpw_r,0,0,nprow)
                         iccol = indxg2p(iw,numpw_c,0,0,npcol)
                         iproc = icrow*npcol+iccol
                         if(myrow==icrow .and. mycol==iccol) then
                            ilrow=indxg2l(jw,numpw_r,0,0,nprow)
                            ilcol=indxg2l(iw,numpw_c,0,0,npcol)
                            rsca=vmat_states(ilrow,ilcol)
                         endif
                         call mp_barrier!ATTENZIONE
                         call mp_bcast(rsca,iproc)
                         tmpspacei(1:ngm_max,iw-(iiw-1)*n_set)=tmpspacei(1:ngm_max,iw-(iiw-1)*n_set)+&
 & rsca*tmpspacej(1:ngm_max,jw-(jjw-1)*n_set)
#endif
                      endif
                   enddo
                endif
             enddo
          enddo
          do iw=(iiw-1)*n_set+1,min(iiw*n_set,numpw)
             if(on_cor(iw) /= 0) then
                CALL davcio(tmpspacei(:,iw-(iiw-1)*n_set), max_ngm*2,iungprod_on,on_cor(iw),1)
             endif
          enddo
       enddo




       close(iungprod)
       close(iungprod_on)


!for compatibility write orthonorm file


       if(ionode) then
          allocate(tmp_str(numpw_on))
          iunu = find_free_unit()

          open(unit=iunu,file=trim(prefix)//'.orthonorm',status='unknown',form='unformatted')
          write(iunu) numpw_on

          do  i=1,numpw_on
             tmp_str(:)=0.d0
             tmp_str(i)=1.d0
             write(iunu) tmp_str
          enddo

          close(iunu)
          deallocate(tmp_str)
       endif

       write(stdout,*) 'TOTAL NUMBER OF ORTHONORMAL PRODUCTS:', numpw_on
       call flush_unit(stdout)



!calculate terms <\tilde{w}_j|w_i>

       write(stdout,*) 'ATTENZIONE1'
       call flush_unit(stdout)


!update overlap matrix omat
       deallocate(tmpspacei,tmpspacej)

        write(stdout,*) 'ATTENZIONE2'
       call flush_unit(stdout)


       if(.not.l_pmatrix) then
          allocate(omat_tmp(numpw_vvc,numpw))
       else
          numpw_vvc_r=n_r
          numpw_vvc_c=n_c
          numpw_vvc_dimr=newdimr
          numpw_vvc_dimc=newdimc
          allocate(omat_tmp(numpw_vvc_dimr,numpw_dimc))
       endif


       if(.not.l_pmatrix) then
          omat_tmp(1:numpw_vvc, 1:numpw)=omat(1:numpw_vvc,1:numpw)
       else
#ifdef __SCALAPACK
          write(stdout,*) 'ATTENZIONE3',numpw_vvc,numpw,numpw_vvc_r,numpw_vvc_dimr,numpw_c,numpw_dimc,n_r,n_c
          call flush_unit(stdout)


          do jw=1,numpw_vvc
             do iw=1,numpw
                icrow = indxg2p(jw,n_r,0,0,nprow)
                iccol = indxg2p(iw,n_c,0,0,npcol)
                iproc=icrow*npcol+iccol
                if(myrow==icrow .and. mycol==iccol) then
                   ilrow=indxg2l(jw,n_r,0,0,nprow)
                   ilcol=indxg2l(iw,n_c,0,0,npcol)
                   rsca=omat(ilrow,ilcol)
                endif
                call mp_barrier!ATTENZIONE
                call mp_bcast(rsca,iproc)

                icrow = indxg2p(jw,numpw_vvc_r,0,0,nprow)
                iccol = indxg2p(iw,numpw_c,0,0,npcol)
                if(myrow==icrow .and. mycol==iccol) then
                   ilrow=indxg2l(jw,numpw_vvc_r,0,0,nprow)
                   ilcol=indxg2l(iw,numpw_c,0,0,npcol)
                   omat_tmp(ilrow,ilcol)=rsca
                endif
             enddo
          enddo
#endif
       endif

        write(stdout,*) 'ATTENZIONE4'
       call flush_unit(stdout)


       if(.not. l_pmatrix) then
          do iw=1,numpw
             if(on_cor(iw)/=0) then
                vmat(:,on_cor(iw))=vmat(:,iw)
             endif
          enddo
     else
#ifdef __SCALAPACK

        do iw=1,numpw
           if(on_cor(iw) /=0) then

              do jw=1,numpw
                 icrow = indxg2p(jw,numpw_r,0,0,nprow)
                 iccol = indxg2p(iw,numpw_c,0,0,npcol)
                 iproc = icrow*npcol+iccol
                 if(myrow==icrow .and. mycol==iccol) then
                    ilrow=indxg2l(jw,numpw_r,0,0,nprow)
                    ilcol=indxg2l(iw,numpw_c,0,0,npcol)
                    rsca=vmat_states(ilrow,ilcol)
                 endif


                 call mp_barrier!ATTENZIONE
                 call mp_bcast(rsca,iproc)



                 icrow = indxg2p(jw,numpw_r,0,0,nprow)
                 iccol = indxg2p(on_cor(iw),numpw_c,0,0,npcol)
                 if(myrow==icrow .and. mycol==iccol) then
                    ilrow=indxg2l(jw,numpw_r,0,0,nprow)
                    ilcol=indxg2l(on_cor(iw),numpw_c,0,0,npcol)
                    vmat(ilrow,ilcol)=rsca
                 endif
              enddo
           endif
        enddo

#endif
     endif

      write(stdout,*) 'ATTENZIONE5'
       call flush_unit(stdout)


     if(.not.l_pmatrix) then
        call dgemm('N','N',numpw_vvc,numpw_on, numpw,1.d0,omat_tmp,numpw_vvc,vmat,numpw,0.d0,omat,numpw_vvc)
     else
#ifdef __SCALAPACK

!pdgemm
! A = omat_tmp
             desc_a(1)=1
             desc_a(2)=icontxt
             desc_a(3)=numpw_vvc
             desc_a(4)=numpw
             desc_a(5)=numpw_vvc_r
             desc_a(6)=numpw_c
             desc_a(7)=0
             desc_a(8)=0
             desc_a(9)=numpw_vvc_dimr

! B = vmat
             desc_b(1)=1
             desc_b(2)=icontxt
             desc_b(3)=numpw
             desc_b(4)=numpw
             desc_b(5)=numpw_r
             desc_b(6)=n_c
             desc_b(7)=0
             desc_b(8)=0
             desc_b(9)=numpw_dimr


! C = omat
             desc_c(1)=1
             desc_c(2)=icontxt
             desc_c(3)=numpw_vvc
             desc_c(4)=numpw_on
             desc_c(5)=n_r
             desc_c(6)=n_c
             desc_c(7)=0
             desc_c(8)=0
             desc_c(9)=newdimr

             call pdgemm('N','N',numpw_vvc,numpw_on, numpw,1.d0,omat_tmp,1,1,desc_a,vmat,1,1,desc_b,0.d0,omat,&
                  &1,1,desc_c)

             write(stdout,*) 'ATTENZIONE6'
             call flush_unit(stdout)


#endif
     endif
     deallocate(omat_tmp)

!update dimension of reduced orthonormal products of wanniers
     numpw_original=numpw
     numpw=numpw_on



     deallocate(eigen, on_cor)
     deallocate(vmat_states)



     return
   end subroutine coulomb_analysis
