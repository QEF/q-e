! FOR GWW
!
! Author: P. Umari
!
 subroutine wannier_pmat_terms_ggrid(numw_prod_dimr,numw_prod_dimc,omat,n_set,ecutoff)

!this subroutine
!calculates the products <wiwj(r_1)|exp(iGX) | w'ij'(r_2)>
!in  wafefunctions G space
!only for the overlapping products of wannier wfcs
!it calculates also the overlaps <wiwj|wi'wj'>
!now w_prod is set up in routine distance_products
!it also writes the elements \tilde{w^P_i} (G=0) on file
!#ifdef __GWW

  USE io_global,            ONLY : stdout, ionode, ionode_id
  USE io_files,             ONLY : find_free_unit, prefix, diropn
  use mp_global,            ONLY : nproc_pool, me_pool
  USE mp,                   ONLY : mp_bcast, mp_sum
  USE kinds,                ONLY : DP
  USE gvect
  USE fft_base,             ONLY : dfftp
  USE basis
  USE klist
  USE constants,            ONLY : e2, pi, tpi, fpi
  USE cell_base,            ONLY : at, alat, tpiba, omega, tpiba2
  USE wvfct,                ONLY : igk, g2kin, npwx, npw, nbnd, nbndx, ecutwfc
  USE wannier_gw



  implicit none

!  COMPLEX(kind=DP), INTENT(inout) :: pmat(3,numw_prod,numw_prod)
  INTEGER, INTENT(in) :: numw_prod_dimr,numw_prod_dimc
  REAL(kind=DP), INTENT(inout)  :: omat(numw_prod_dimr,numw_prod_dimc)
!  INTEGER, INTENT(in) :: gtable(npwx,3)!G-->G+1 correspondence
  INTEGER, INTENT(in) :: n_set  !defines the number of states to be read from disk at the same time
  REAL(kind=DP), INTENT(in) :: ecutoff!cutoff in Rydberg for g sum


  !  --- Internal definitions ---

  INTEGER :: iw,jw,iiw,jjw,jw_begin
  INTEGER :: i,j,ig
  INTEGER :: iungprod
  COMPLEX(kind=DP), ALLOCATABLE :: tmpgi(:,:),tmpgj(:,:)
  LOGICAL :: match
  INTEGER :: ndist
  INTEGER :: ndistance!this is a function
  INTEGER :: ix,iy,iz
  INTEGER :: iunprod
  TYPE(wannier_product) :: wiwj
  INTEGER :: rspacel(3)
  INTEGER :: igk0(npwx)
  REAL(kind=dp) :: g2kin_bp(npwx)
  INTEGER :: npw0
  LOGICAL :: exst
  INTEGER :: mdir

  INTEGER :: iun_g
  REAL(kind=DP), ALLOCATABLE :: wpg(:)
  REAL(kind=DP) :: sca
  INTEGER :: ngm_max
  REAL(kind=DP), ALLOCATABLE :: omat_tmp(:,:)
  INTEGER :: iw_min,iw_max,jw_min,jw_max
  INTEGER :: numw_prod_r,numw_prodc
  INTEGER :: icrow,iccol,iproc,ilrow,ilcol
#ifdef __SCALAPACK
  INTEGER, EXTERNAL :: indxg2p,indxg2l
#endif
  write(stdout,*) 'ATTENZIONE0'
  call flush_unit(stdout)

  if(l_pmatrix) then
     numw_prod_r=ceiling(real(numw_prod)/real(max(nprow,npcol)))
     numw_prodc=ceiling(real(numw_prod)/real(max(nprow,npcol)))
  endif

!determine ngm_max
   ngm_max=0
   do ig=1,ngm
      if(gg(ig)*tpiba2 >= ecutoff) exit
      ngm_max=ngm_max+1
   enddo

   write(stdout,*) 'NGM MAX:', ngm_max, ngm
   call flush_unit(stdout)


  rspacel(1)=dfftp%nr1
  rspacel(2)=dfftp%nr2
  rspacel(3)=dfftp%nr3

!  if(.not.lnonorthogonal) pmat(1:3,1:numw_prod,1:numw_prod)=(0.d0,0.d0)
  omat(:,:)=0.d0

  write(stdout,*) 'Routine: wannier_pmat_terms',n_set
 call flush_unit(stdout)

  allocate(wpg(numw_prod))
  if(ionode) then
      iun_g =  find_free_unit()
      open( unit= iun_g, file=trim(prefix)//'.gzero', status='unknown',form='unformatted')
      write(iun_g) numw_prod
   endif


   CALL gk_sort(xk(1,1),ngm,g,ecutwfc/tpiba2, &
              &    npw0,igk0,g2kin_bp)
   iungprod = find_free_unit()


   allocate(tmpgi(max_ngm,n_set),tmpgj(max_ngm,n_set))
   allocate(omat_tmp(n_set,n_set))

   CALL diropn( iungprod, 'wiwjwfc', max_ngm*2, exst )

   write(stdout,*) 'ATTENZIONE1'
   call flush_unit(stdout)

!allocate and reads wannier product descriptors
  ALLOCATE(w_prod(numw_prod))
  if(ionode) then
     iunprod = find_free_unit()
     open( unit= iunprod, file=trim(prefix)//'.wiwjprod', status='unknown',form='unformatted')
     rewind(iunprod)
     do iw=1,numw_prod
        read(iunprod)  wiwj
        w_prod(iw)=wiwj
     enddo
     close(iunprod)
  endif
  write(stdout,*) 'ATTENZIONE2'
   call flush_unit(stdout)

#ifdef __MPI
  do iw=1,numw_prod
     call mp_bcast(w_prod(iw)%i,ionode_id)
     call mp_bcast(w_prod(iw)%j,ionode_id)
     call mp_bcast(w_prod(iw)%center(:),ionode_id)
     call mp_bcast(w_prod(iw)%radius(:),ionode_id)
  enddo
#endif
  do iiw=1,numw_prod/n_set+1
     write(stdout,*) 'IIW', iiw
     call flush_unit(stdout)
      do iw=(iiw-1)*n_set+1,min(iiw*n_set,numw_prod)

         CALL davcio(tmpgi(:,iw-(iiw-1)*n_set),max_ngm*2,iungprod,iw,-1)

         sca=0.d0
         if(gstart==2) sca=dble(tmpgi(1,iw-(iiw-1)*n_set))
         call mp_sum(sca)
         if(ionode) write(iun_g) sca


      enddo
      write(stdout,*) 'ATTENZIONE3'
      call flush_unit(stdout)

      iw_min=(iiw-1)*n_set+1
      iw_max=min(iiw*n_set,numw_prod)

      do jjw=iiw,numw_prod/n_set+1
         write(stdout,*) 'JJW', jjw
         call flush_unit(stdout)

         do jw=(jjw-1)*n_set+1,min(jjw*n_set,numw_prod)
            CALL davcio(tmpgj(:,jw-(jjw-1)*n_set),max_ngm*2,iungprod,jw,-1)

         enddo
         jw_min=(jjw-1)*n_set+1
         jw_max=min(jjw*n_set,numw_prod)
!!!!!!!!!!!!!!!!!!!!!!!!!!!
!uses  blas routine


         if(gstart==2) tmpgj(1,:)=0.5d0*tmpgj(1,:)
!         call zgemm('C','N',n_set,n_set,ngm_max,(1.d0,0.d0),tmpgi,max_ngm,tmpgj,max_ngm,(0.d0,0.d0),omat_tmp,n_set)
!         call mp_sum(omat_tmp(:,:))
!         do iw=iw_min,iw_max
!            do jw=jw_min,jw_max
!               omat(iw,jw)=2.d0*dble(omat_tmp(iw-iw_min+1,jw-jw_min+1))
!               omat(jw,iw)=omat(iw,jw)
!            enddo
!         enddo
         call dgemm('T','N',n_set,n_set,2*ngm_max,2.d0,tmpgi,2*max_ngm,tmpgj,2*max_ngm,0.d0,omat_tmp,n_set)
         call mp_sum(omat_tmp(:,:))
          do iw=iw_min,iw_max
            do jw=jw_min,jw_max
               if(.not.l_pmatrix) then
                  omat(iw,jw)=omat_tmp(iw-iw_min+1,jw-jw_min+1)
                  omat(jw,iw)=omat(iw,jw)
               else
#ifdef __SCALAPACK
                  icrow = indxg2p(jw,numw_prod_r,0,0,nprow)
                  iccol = indxg2p(iw,numw_prodc,0,0,npcol)
                  if(myrow==icrow .and. mycol==iccol) then
                     ilrow=indxg2l(jw,numw_prod_r,0,0,nprow)
                     ilcol=indxg2l(iw,numw_prodc,0,0,npcol)
                     omat(ilrow,ilcol)=omat_tmp(iw-iw_min+1,jw-jw_min+1)
                  endif
                  icrow = indxg2p(iw,numw_prod_r,0,0,nprow)
                  iccol = indxg2p(jw,numw_prodc,0,0,npcol)
                  if(myrow==icrow .and. mycol==iccol) then
                     ilrow=indxg2l(iw,numw_prod_r,0,0,nprow)
                     ilcol=indxg2l(jw,numw_prodc,0,0,npcol)
                     omat(ilrow,ilcol)=omat_tmp(iw-iw_min+1,jw-jw_min+1)
                  endif
#endif
               endif
            enddo
         enddo
      enddo
   enddo
  DEALLOCATE(tmpgi,tmpgj)
  deallocate(omat_tmp)
  CLOSE(iungprod)

  if(ionode) close(iun_g)
  deallocate(wpg)

!#endif
  return
end subroutine wannier_pmat_terms_ggrid
