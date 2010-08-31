! FOR GWW
!
! Author: P. Umari
!
  subroutine  set_wannier_P(numpw,numpw_vvc,w_prod, qmat,dimr,dimc,tresh,w_P,wp_dimr, wp_r,n_r,n_c)

! this subroutine finds the non zero overlap (>tresh) of
! orthonormalized and localized products of wanniers w_P
! with wanniers w_i*w_j

   USE kinds, ONLY : DP
   USE wannier_gw, ONLY : wannier_product, wannier_P,  l_pmatrix, npcol, nprow,icontxt,myrow,mycol
   USE io_global, ONLY : stdout,ionode
   USE mp, ONLY : mp_bcast
   USE mp_global, ONLY : nproc, mpime

   implicit none

   INTEGER :: numpw!number of products of wannier
   INTEGER :: numpw_vvc!total number of direct products of wannier
   TYPE(wannier_product) :: w_prod(numpw_vvc)!descriptors of NOT orthonormal products of wannier
   INTEGER, INTENT(in) :: dimr,dimc!dimensions of qmat
   REAL(kind=DP) :: qmat(dimr,dimc)!overlaps <w_P|w_i*w_j>
   REAL(kind=DP) :: tresh!treshold for absolute value of products
   TYPE(wannier_P) :: w_P(numpw)!descriptors of  orthonormal products of wannier
   INTEGER, INTENT(in)  :: wp_dimr!dimension of w_p
   INTEGER, INTENT(in)  :: wp_r !multiplicity of w_p on nproc for parallel execution
   INTEGER, INTENT(in) :: n_r,n_c!in parallel case periodicity on rows and columns

!internal variables
   INTEGER :: iw,jw,ip,ii
   INTEGER, ALLOCATABLE :: ij(:,:)
   REAL(kind=DP), ALLOCATABLE :: o(:)
   INTEGER :: noverlap
   INTEGER :: icrow,iccol,ilrow,ilcol,iproc
   REAL(kind=DP) :: rsca
#ifdef __SCALAPACK
   INTEGER, EXTERNAL :: indxg2p,indxg2l
#endif
   REAL(kind=DP), ALLOCATABLE :: sndbuf(:,:)
   INTEGER, ALLOCATABLE :: ipos(:), ibuf(:,:)
   REAL(kind=DP), ALLOCATABLE :: rvect(:)

   write(stdout,*) 'set_wannier_P',dimr,dimc,n_r,n_c!ATTENZIONE
   write(stdout,*) 'NUMPW NUMPW_VCC', numpw, numpw_vvc
   allocate(ij(2,numpw_vvc),o(numpw_vvc))

   if(.not.l_pmatrix) then
!ONLY IONODE_ID is allocating this array
      if(ionode) then
         do iw=1,numpw!loop on w_P
            noverlap=0
            do jw=1,numpw_vvc!loop on w_prod
               rsca=qmat(jw,iw)
               if(abs(rsca) >= tresh) then
                  noverlap=noverlap+1
                  ij(1,noverlap)=w_prod(jw)%i
                  ij(2,noverlap)=w_prod(jw)%j
                  o(noverlap)=rsca
! write(*,*) 'w_P', iw,ij(1,noverlap),ij(2,noverlap),noverlap,qmat(jw,iw)
               else
                  write(stdout,*) 'Set_wannier_P eliminated:', iw, jw,rsca
               endif
            enddo
            w_P(iw)%numij=noverlap
            allocate(w_P(iw)%ij(2,noverlap))
            allocate(w_P(iw)%o(noverlap))
            w_P(iw)%ij(1,1:noverlap)=ij(1,1:noverlap)
            w_P(iw)%ij(2,1:noverlap)=ij(2,1:noverlap)
            w_P(iw)%o(1:noverlap)=o(1:noverlap)
         enddo
      endif
   else
#ifdef __SCALAPACK
      do iw=1,numpw!loop on w_P
         noverlap=0
         allocate(sndbuf(numpw_vvc,nprow*npcol),ibuf(numpw_vvc,nprow*npcol))
         allocate(ipos(nprow*npcol),rvect(numpw_vvc))
         ipos(:)=0
         do jw=1,numpw_vvc
            icrow = indxg2p(jw,n_r,0,0,nprow)
            iccol = indxg2p(iw,n_c,0,0,npcol)
            iproc=icrow*npcol+iccol
            ipos(iproc+1)=ipos(iproc+1)+1
            ibuf(ipos(iproc+1),iproc+1)=jw
            if(myrow==icrow .and. mycol==iccol) then
               ilrow=indxg2l(jw,n_r,0,0,nprow)
               ilcol=indxg2l(iw,n_c,0,0,npcol)
               sndbuf(ipos(iproc+1),iproc+1)=qmat(ilrow,ilcol)
            endif
         enddo
         do ip=1,nprow*npcol
            if(ipos(ip)>0) then
               call mp_bcast(sndbuf(1:ipos(ip),ip),ip-1)
               do ii=1,ipos(ip)
                  jw=ibuf(ii,ip)
                  rvect(jw)=sndbuf(ii,ip)
               enddo
            endif
         enddo
         do jw=1,numpw_vvc
            if(abs(rvect(jw)) >= tresh) then
               noverlap=noverlap+1
               ij(1,noverlap)=w_prod(jw)%i
               ij(2,noverlap)=w_prod(jw)%j
               o(noverlap)=rvect(jw)
            else
               write(stdout,*) 'Set_wannier_P eliminated:', iw, jw,rvect(jw)
            endif
         enddo
         if(indxg2p(iw,wp_r,0,0,nproc)==mpime) then
            ilrow=indxg2l(iw,wp_r,0,0,nproc)
            w_P(ilrow)%numij=noverlap
            allocate(w_P(ilrow)%ij(2,noverlap))
            allocate(w_P(ilrow)%o(noverlap))
            w_P(ilrow)%ij(1,1:noverlap)=ij(1,1:noverlap)
            w_P(ilrow)%ij(2,1:noverlap)=ij(2,1:noverlap)
            w_P(ilrow)%o(1:noverlap)=o(1:noverlap)
         endif
         deallocate(sndbuf,ibuf,ipos,rvect)
      enddo
#endif
   endif

   return
 end subroutine set_wannier_P
