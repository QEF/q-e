
 subroutine allocate_bp_efield
!this subroutine allocate memory for the Berry's phase electric field 
   USE bp
   USE ions_base,  ONLY : nat
   USE gvect, ONLY : ngm_g

   implicit none
   
   if(lberry.or.lelfield) then
      allocate(mapgp_global(ngm_g,3))
      allocate(mapgm_global(ngm_g,3))
      allocate(forces_bp_efield(3,nat))
   endif

   l_el_pol_old=.false.
   el_pol_acc=0.d0

   return
 end subroutine allocate_bp_efield


 subroutine deallocate_bp_efield
!this subroutine allocate memory for the Berry's phase electric field
   USE bp

   implicit none

   if(lberry.or.lelfield) then
      deallocate(mapgp_global)
      deallocate(mapgm_global)
      deallocate(forces_bp_efield)
      if(allocated(nx_el)) deallocate(nx_el)
   endif

   return
 end subroutine deallocate_bp_efield

 
 subroutine bp_global_map
!this subroutine sets up the global correspondence map G+1 and G-1

    USE kinds,                ONLY : DP
    USE io_global,            ONLY : stdout
    USE mp,                   ONLY : mp_sum
    USE bp
    USE gvect,                ONLY : ngm_g, g, ngm
    USE fft_base,             ONLY : dfftp
    USE cell_base,            ONLY : at
    USE gvect,   ONLY : ig_l2g

    implicit none

    INTEGER :: ig, mk1,mk2,mk3, idir, imk(3)
    INTEGER, ALLOCATABLE :: ln_g(:,:,:)
    INTEGER, ALLOCATABLE :: g_ln(:,:)

!set up correspondence ln_g ix,iy,iz ---> global g index in
! (for now...) coarse grid
!and inverse realtion global g (coarse) to ix,iy,iz

    allocate(ln_g(-dfftp%nr1:dfftp%nr1,-dfftp%nr2:dfftp%nr2,-dfftp%nr3:dfftp%nr3))
    allocate(g_ln(3,ngm_g))
    
    ln_g(:,:,:)=0!it means also not found
    do ig=1,ngm
       mk1=nint(g(1,ig)*at(1,1)+g(2,ig)*at(2,1)+g(3,ig)*at(3,1))
       mk2=nint(g(1,ig)*at(1,2)+g(2,ig)*at(2,2)+g(3,ig)*at(3,2))
       mk3=nint(g(1,ig)*at(1,3)+g(2,ig)*at(2,3)+g(3,ig)*at(3,3))
       ln_g(mk1,mk2,mk3)=ig_l2g(ig)
    enddo
    call mp_sum(ln_g(:,:,:))
    

    g_ln(:,:)= 0!it means also not found
    do ig=1,ngm
       mk1=nint(g(1,ig)*at(1,1)+g(2,ig)*at(2,1)+g(3,ig)*at(3,1))
       mk2=nint(g(1,ig)*at(1,2)+g(2,ig)*at(2,2)+g(3,ig)*at(3,2))
       mk3=nint(g(1,ig)*at(1,3)+g(2,ig)*at(2,3)+g(3,ig)*at(3,3))
       g_ln(1,ig_l2g(ig))=mk1
       g_ln(2,ig_l2g(ig))=mk2
       g_ln(3,ig_l2g(ig))=mk3
    enddo
    call mp_sum(g_ln(:,:))
    
!loop on direction
    do idir=1,3
!for every g on global array find G+1 and G-1 and put on
       do ig=1,ngm_g
          imk(:)=g_ln(:,ig)

          
          imk(idir)=imk(idir)+1
!table array
          mapgp_global(ig,idir)=ln_g(imk(1),imk(2),imk(3))
          imk(idir)=imk(idir)-2
          mapgm_global(ig,idir)=ln_g(imk(1),imk(2),imk(3))
       enddo
    enddo
    deallocate(ln_g,g_ln)


    return

  end subroutine bp_global_map
