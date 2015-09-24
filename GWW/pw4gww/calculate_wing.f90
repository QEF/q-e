!
! Copyright (C) 2001-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!

subroutine calculate_wing(n_set, orthonorm)

!this subroutine calculate the terms
!\Sum_G <G|\tilde{w^P_i}>\epsilon(G'=0, G, iw)
! it requires the file .e_head 
  USE io_global,            ONLY : stdout, ionode, ionode_id
  USE io_files,             ONLY : prefix, tmp_dir, diropn
  USE kinds,                ONLY : DP
  USE wannier_gw
  USE mp,                   ONLY : mp_bcast, mp_sum
  USE gvect,                ONLY : mill, ngm, gstart,g,ngm_g, ig_l2g
  USE cell_base,            ONLY : tpiba
  USE mp_wave, ONLY : mergewf,splitwf
  USE mp_global, ONLY : intra_pool_comm
  USE mp_world,  ONLY : mpime, nproc, world_comm
  USE wvfct,    ONLY :  npwx, npw
  USE cell_base, ONLY : at,bg

  implicit none

  INTEGER, EXTERNAL :: find_free_unit

  INTEGER, INTENT(in)  :: n_set  !defines the number of states to be read from disk at the same time
  INTEGER, INTENT(in)  :: orthonorm!if ==1 opens orthonormalized products of wannier file, if ==2 reduced one


  INTEGER iun, iungprod
  INTEGER :: n_g, ngm_k
  INTEGER :: ig, igg, iw, iiw, i, ii, ipol
  REAL(kind=DP) :: omega_g
  REAL(kind=DP), ALLOCATABLE :: freqs(:)
  COMPLEX(kind=DP), ALLOCATABLE :: e_head(:,:,:), e_head_g0(:)
  LOGICAL :: exst
  COMPLEX(kind=DP), ALLOCATABLE :: tmpspacei(:,:)
  REAL(kind=DP), ALLOCATABLE :: wing(:,:,:), wing_c(:,:,:)
  REAL(kind=DP) :: sca
  REAL(kind=DP), ALLOCATABLE :: fact(:)
  REAL(kind=DP) :: qq
  INTEGER :: npwx_g
  INTEGER, ALLOCATABLE :: k2g_ig_l2g(:)

  npwx_g=npwx
  call mp_sum(npwx_g,world_comm)

!read file .e_head

  write(stdout,*) 'Routine calculate_wing'
  FLUSH(stdout)
  allocate(fact(ngm))
  allocate(k2g_ig_l2g(ngm))

  fact(:)=0.d0
  if(gstart==2) fact(1)=0.d0
  do ig=gstart,npw
!     qq = g(1,ig)**2.d0 + g(2,ig)**2.d0 + g(3,ig)**2.d0
!     fact(ig)=1.d0/tpiba/dsqrt(qq)
     fact(ig)=dsqrt(vg_q(ig))
  end do


  call ktogamma_ig_l2g ( k2g_ig_l2g, at, bg )

  write(stdout,*) 'ATT0.1'
  FLUSH(stdout)

  
  
  if(ionode) then
      iun =  find_free_unit()
      open( unit= iun, file=trim(tmp_dir)//'/_ph0/'//trim(prefix)//'.e_head', status='old',form='unformatted')
      read(iun) n_g
      read(iun) omega_g
   endif


   call mp_bcast(n_g, ionode_id,world_comm)
   call mp_bcast(omega_g, ionode_id,world_comm)
   allocate(freqs(n_g+1))

   if(ionode) then
      read(iun) freqs(1:n_g+1)
      read(iun) ngm_k
   endif

  write(stdout,*) 'ATT0.2'
  FLUSH(stdout)


   call mp_bcast(freqs(:), ionode_id,world_comm)
   call mp_bcast(ngm_k, ionode_id,world_comm)

   allocate(e_head_g0(ngm_k))
   allocate(e_head(npw, n_g+1,3)) 
   e_head(:,:,:) = (0.d0,0.d0)  
   do ii=1,n_g+1
      do ipol=1,3
         e_head_g0(:)=(0.d0,0.d0)
         if(ionode) read(iun) e_head_g0(1:ngm_k)
         call splitwf(e_head(:, ii,ipol),e_head_g0,npw,k2g_ig_l2g,mpime,nproc,ionode_id,intra_pool_comm) 
      enddo
   enddo
    if(ionode) close(iun)

   write(stdout,*) 'ATT1'
   FLUSH(stdout)

   deallocate(e_head_g0)
!loop on n_set groups
   write(stdout,*) 'ATT2'
   FLUSH(stdout)


   allocate(tmpspacei(max_ngm,n_set))
   iungprod = find_free_unit()

   if(orthonorm==0) then
      CALL diropn( iungprod, 'wiwjwfc', max_ngm*2, exst )
   else if(orthonorm==1) then
      CALL diropn( iungprod, 'wiwjwfc_on', max_ngm*2, exst )
   else
      CALL diropn( iungprod, 'wiwjwfc_red', max_ngm*2, exst )
   endif
   allocate(wing(numw_prod, n_g+1,3))
   wing(:,:,:)=0.d0
   !allocate(wing_c(numw_prod, n_g+1,3))
   !wing_c(:,:,:)=0.d0



   do iiw=1,ceiling(real(numw_prod)/real(n_set))
!read states
      do iw=(iiw-1)*n_set+1,min(iiw*n_set,numw_prod)
         CALL davcio(tmpspacei(:,iw-(iiw-1)*n_set),max_ngm*2,iungprod,iw,-1)
      enddo


      write(stdout,*) 'ATT3'
      FLUSH(stdout)


 !loop on states
      do iw=(iiw-1)*n_set+1,min(iiw*n_set,numw_prod)
         do ipol=1,3
            do i=1, n_g+1
               sca=0.d0
               do ig=1,max_ngm
!                  sca=sca+2.d0*real(tmpspacei(ig,iw-(iiw-1)*n_set)*conjg(e_head(ig,i,ipol)))*fact(ig)
                   sca=sca+2.d0*dble((tmpspacei(ig,iw-(iiw-1)*n_set))*conjg(e_head(ig,i,ipol)))!*fact(ig)!ATTENZIONE
               enddo
               call mp_sum(sca,world_comm)
               wing(iw,i,ipol)=sca


            enddo
         enddo
         
      enddo
   enddo

    write(stdout,*) 'ATT4'
   FLUSH(stdout)
   

!write terms on file

    if(ionode) then
       iun =  find_free_unit()
       open( unit= iun, file=trim(tmp_dir)//trim(prefix)//'.wing', status='unknown',form='unformatted')
       write(iun) n_g
       write(iun) omega_g
       write(iun) numw_prod
       do ipol=1,3
          do i=1,n_g+1
             write(iun) wing(1:numw_prod,i,ipol)
          enddo
!          do i=1,n_g+1
!             write(iun) wing_c(1:numw_prod,i,ipol)
!          enddo
       enddo
       close(iun)
   endif

   deallocate(tmpspacei)
   close(iungprod)
   deallocate(fact)
   deallocate(freqs)
  ! if(ionode) deallocate (e_head_g)
   deallocate(e_head)
   deallocate(wing)
   deallocate(k2g_ig_l2g)
  return
end subroutine calculate_wing


   !-----------------------------------------------------------------------
   SUBROUTINE ktogamma_ig_l2g ( k2g_ig_l2g, at, bg )
   !----------------------------------------------------------------------
   !
   !     This routine generates all the reciprocal lattice vectors
   !     contained in the sphere of radius gcutm. Furthermore it
   !     computes the indices nl which give the correspondence
   !     between the fft mesh points and the array of g vectors.
   !
   USE gvect,              ONLY : ig_l2g, g, gg, ngm, ngm_g, gcutm, &
                                  mill,  nl, gstart
   USE fft_base,           ONLY : dfftp, dffts
!                                                                                                                                                           
   USE kinds,              ONLY : DP
   USE constants,          ONLY : eps8





   IMPLICIT NONE
   !
   INTEGER, INTENT(out) :: k2g_ig_l2g(ngm)
   REAL(DP), INTENT(IN) ::  at(3,3), bg(3,3)
   !     here a few local variables
   !
   REAL(DP) ::  t (3), tt
   INTEGER :: ngm_, n1, n2, n3, n1s, n2s, n3s
   !
   REAL(DP), ALLOCATABLE :: g2sort_g(:)
   ! array containing all g vectors, on all processors: replicated data
   INTEGER, ALLOCATABLE :: mill_g(:,:), mill_unsorted(:,:)
   ! array containing all g vectors generators, on all processors:
   !     replicated data
   INTEGER, ALLOCATABLE :: igsrt(:)
   !
   INTEGER :: ni, nj, nk, i, j, k, ipol, ng, igl, indsw,ig
   !
   ! counters
   !
   !    set the total number of fft mesh points and and initial value of gg
   !    The choice of gcutm is due to the fact that we have to order the
   !    vectors after computing them.
   !
   !
   !     set d vector for unique ordering
   !
   !    and computes all the g vectors inside a sphere
   !
   ALLOCATE( mill_g( 3, ngm_g*3 ),mill_unsorted( 3, ngm_g*3 ) )
   ALLOCATE( igsrt( ngm_g*3 ) )
   ALLOCATE( g2sort_g( ngm_g*3 ) )
   g2sort_g(:) = 1.0d20
   !
   ! save present value of ngm 
   !
  
   !
   ngm_ = 0
  
   !
   ! max miller indices (same convention as in module stick_set)
   !
   ni = (dfftp%nr1-1)/2
   nj = (dfftp%nr2-1)/2
   nk = (dfftp%nr3-1)/2
   !
   iloop: DO i = -ni, ni
      !
    
      jloop: DO j = -nj, nj
         !
    
         kloop: DO k = -nk, nk
                !
    
            t(:) = i * bg (:,1) + j * bg (:,2) + k * bg (:,3)
            tt = sum(t(:)**2)
            IF (tt <= gcutm) THEN
               ngm_ = ngm_ + 1
               mill_unsorted( :, ngm_ ) = (/ i,j,k /)
               IF ( tt > eps8 ) THEN
                  g2sort_g(ngm_) = tt
               ELSE
                  g2sort_g(ngm_) = 0.d0
               ENDIF
            ENDIF
         ENDDO kloop
      ENDDO jloop
   ENDDO iloop


   igsrt(1) = 0
   CALL hpsort_eps( ngm_, g2sort_g, igsrt, eps8 )
   mill_g(1,1:ngm_) = mill_unsorted(1,igsrt(1:ngm_))
   mill_g(2,1:ngm_) = mill_unsorted(2,igsrt(1:ngm_))
   mill_g(3,1:ngm_) = mill_unsorted(3,igsrt(1:ngm_))
   DEALLOCATE( g2sort_g, igsrt, mill_unsorted )



   do ig=1,ngm
      do ng=1,ngm_
         if(mill_g(1,ng)==mill(1,ig).and.mill_g(2,ng)==mill(2,ig).and.mill_g(3,ng)==mill(3,ig)) then
            k2g_ig_l2g(ig)=ng
            exit
         endif
      enddo
   enddo

   !
   DEALLOCATE( mill_g )



 END SUBROUTINE ktogamma_ig_l2g
