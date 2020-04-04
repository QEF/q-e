MODULE simple_objects
!this module describes the most important objects

  USE kinds, ONLY : DP

  TYPE bands
!data regarding bands and wavefunctions
     INTEGER :: numv !number of valence states (those considered for excitons only)
     INTEGER :: numc !number of conduction states
     INTEGER :: num!numv+numc
     INTEGER :: ntot_e!dimension of global to all k, basis for KS states
     INTEGER :: nk!total number of k, points
     INTEGER :: nk_loc!local number of k points
     INTEGER :: ik_first!first local k point
     INTEGER :: ik_last!last local k point
     REAL(kind=DP) :: scissor!scissor in Ry a.u.
     REAL(kind=DP) :: ene(4)!tot,diagonal,exchange,direct
     REAL(kind=DP), DIMENSION(:,:), POINTER :: k!coordinates of local k points in atomic units 
     COMPLEX(kind=DP), DIMENSION(:,:,:), POINTER :: omat!overlap matrix (ntot_e,num,nk_loc)
     REAL(kind=DP), DIMENSION(:,:), POINTER :: en_v!KS or GW  energies of valence states (numv,nk_loc) 
     REAL(kind=DP), DIMENSION(:,:), POINTER :: en_c!KS or GW  energies of valence states (numc,nk_loc)

  END TYPE bands

  TYPE exc
!this describes an electron-hole excitons
      INTEGER :: numv !number of valence states (those considered for excitons only)
      INTEGER :: numc !number of conduction states
      INTEGER :: num!numv+numc           
      INTEGER :: nk!total number of k, points            
      INTEGER :: nk_loc!local number of k points 
      INTEGER :: ik_first!first local k point 
      INTEGER :: ik_last!last local k point
      REAL(kind=DP) :: ene(4)!tot,diagonal,exchange,direct 
      COMPLEX(kind=DP), DIMENSION(:,:,:), POINTER  :: avc!A_vck terms (numv,numc,nk_loc)

  END  TYPE exc

  TYPE product
!this type describes the basis for the products of global wave-function basis vectors
!this is not distributed with k-points parallelization
     INTEGER :: nprod_e!number of product terms
     INTEGER :: ntot_e!dimension of global to all k, basis for KS states
     COMPLEX(kind=DP), DIMENSION(:,:,:), POINTER :: fij!<\mathcal{E_\alpha}|(e_{i}^*e_{j}> terms 
  END TYPE product


  TYPE potential
!this object described a potential matrix in the space of the common basis for products
     INTEGER :: nprod_e!number of product terms 
     COMPLEX(kind=DP), DIMENSION(:,:), POINTER :: vpot
     INTEGER :: nkp(3)!k-points grid equally spaced
     INTEGER :: nk!total number of k, points
     INTEGER, DIMENSION(:,:,:), POINTER :: ijk!correspondence k'-k-->i,j,k (3,k,k')
     COMPLEX(kind=DP), DIMENSION(:,:,:,:,:), POINTER :: vpotq! V_{k'-k} (nprod_e,nprod_e,i,j,k), ijk from 1 to nkp() IMPORTANT
     COMPLEX(kind=DP), DIMENSION(:,:,:,:,:), POINTER :: wpotq! Wc_{k'-k} (nprod_e,nprod_e,i,j,k), ijk from 1 to nkp() IMPORTANT
  END TYPE potential

  TYPE epe
     INTEGER :: ntot_e
     COMPLEX(KIND=DP), DIMENSION(:,:,:), POINTER :: mat
  END TYPE epe
  
  INTERFACE OPERATOR(+)

     MODULE PROCEDURE sum_exc

  END INTERFACE


  INTERFACE OPERATOR(*)

    MODULE PROCEDURE prod_exc
    MODULE PROCEDURE prod_c_exc

  END INTERFACE


  INTERFACE ASSIGNMENT(=)
     MODULE PROCEDURE assign_exc 
  END INTERFACE


  INTERFACE OPERATOR(.hd.)
     MODULE PROCEDURE h_diag
  END INTERFACE
  CONTAINS

!!!!!!!

!      INTEGER :: nk_loc
!      INTEGER :: ik_first
!      INTEGER :: ik_last
!      COMPLEX(kind=DP), DIMENSION(:,:,:), POINTER  :: avc

!  END  TYPE exc

  SUBROUTINE nice_write_exc(bd,simple_in,a,label)
!write in characters an exciton, with label 
    USE input_simple_exc, ONLY : input_options
    USE mp,                   ONLY : mp_sum
      USE mp_world,             ONLY : world_comm
      USE io_files,  ONLY : tmp_dir, prefix
      USE io_global, ONLY : ionode_id, ionode, stdout
      USE mp_world,  ONLY : mpime, nproc

    implicit none
    TYPE(bands) :: bd
    TYPE(input_options) :: simple_in
    TYPE(exc) :: a
    INTEGER :: label

    INTEGER, EXTERNAL :: find_free_unit
    INTEGER :: iun
    CHARACTER(4) :: nfile
    INTEGER :: ik, ikl,iv,ic
    REAL(kind=DP) :: kappa(3)
    COMPLEX(kind=DP), ALLOCATABLE :: avc(:,:)


    allocate(avc(a%numv,a%numc))

    write(nfile,'(4i1)') &
               & label/1000,mod(label,1000)/100,mod(label,100)/10,mod(label,10)


    if(ionode) then
       iun = find_free_unit()
       open( unit=iun, file=trim(tmp_dir)//trim(simple_in%prefix)//'.exciton.'//nfile, status='unknown')
       write(iun,*) a%numv,a%numc,a%nk
       write(iun,*) a%ene(1:4)
    end if
    do ik= 1,a%nk
       kappa=0.d0
       do ikl=a%ik_first,a%ik_last
          if(ikl==ik) then
             kappa(1:3)=bd%k(1:3,ik-a%ik_first+1)
          endif
       enddo
       call mp_sum(kappa,world_comm)
       if(ionode) write(iun,*) kappa(1:3)
    enddo
    do ik= 1,a%nk
       avc=(0.d0,0.d0)
       do ikl=a%ik_first,a%ik_last
          if(ikl==ik) then
             avc(1:a%numv,1:a%numc)=a%avc(1:a%numv,1:a%numc,ik-a%ik_first+1)
          endif
       enddo
       call mp_sum(avc,world_comm)
       if(ionode) then
          do ic=1,a%numc
             do iv=1,a%numv
                write(iun,*) avc(iv,ic)
             enddo
          enddo
       endif
    enddo
    

    if(ionode) close(iun)
    deallocate(avc)
    return
  END SUBROUTINE nice_write_exc
!!!!!!




    subroutine initialize_epe(element)
      implicit none
      TYPE(epe) :: element
      nullify(element%mat)
      return
    end subroutine initialize_epe
    
    subroutine deallocate_epe(element)
      implicit none
      TYPE(epe) :: element
      if(associated(element%mat)) deallocate(element%mat)
      nullify(element%mat)
      return
    end subroutine deallocate_epe

    subroutine read_epe(simple_in,element)

      USE input_simple_exc, ONLY : input_options
      USE mp,                   ONLY : mp_bcast
      USE mp_world,             ONLY : world_comm
      USE io_files,  ONLY : tmp_dir, prefix
      USE io_global, ONLY : ionode_id, ionode, stdout
      USE mp_world,  ONLY : mpime, nproc
      
      implicit none

      TYPE(input_options) :: simple_in
      TYPE(epe) :: element
      INTEGER, EXTERNAL :: find_free_unit
      INTEGER :: iun,i,a
      write(*,*)'epe:opening file'
      if(ionode) then
         iun = find_free_unit()
         open( unit=iun, file=trim(tmp_dir)//trim(simple_in%prefix)//'.epe', status='old',form='unformatted')
         write(*,*)'file opened'
         read(iun) element%ntot_e
      end if
      call mp_bcast(element%ntot_e,ionode_id,world_comm)
      write(*,*)element%ntot_e
      allocate(element%mat(element%ntot_e,element%ntot_e,3))
      if(ionode) then
         do a=1,3
            do i=1,element%ntot_e
               read(iun) element%mat(1:element%ntot_e,i,a)
            end do
         end do
         close(iun)
      end if
      call mp_bcast(element%mat,ionode_id,world_comm)
      write(*,*)'epe:read all'
      return
      
    end subroutine read_epe
    
    subroutine deallocate_bands(bd)
      implicit none
      TYPE(bands) :: bd
      if(associated(bd%k)) deallocate(bd%k)
      if(associated(bd%omat)) deallocate(bd%omat)
      if(associated(bd%en_v)) deallocate(bd%en_v)
      if(associated(bd%en_c)) deallocate(bd%en_c)
      nullify(bd%k)
      nullify(bd%omat)
      nullify(bd%en_v)
      nullify(bd%en_c)
      return
    end subroutine deallocate_bands

    subroutine deallocate_exc(a)
      implicit none
      TYPE(exc) :: a
      if(associated(a%avc)) deallocate(a%avc)
      nullify(a%avc)
      return
    end subroutine deallocate_exc
      
    subroutine deallocate_product(pd)
      implicit none
      TYPE(product) :: pd
      if(associated(pd%fij)) deallocate(pd%fij)
      nullify(pd%fij)
      return
    end subroutine deallocate_product

    subroutine initialize_product(pd)
      implicit none
      TYPE(product) :: pd
      nullify(pd%fij)
    end subroutine initialize_product


    subroutine initialize_potential(pt)
      implicit none
      TYPE(potential) :: pt
      nullify(pt%vpot)
      nullify(pt%vpotq)
      nullify(pt%ijk)
      return
    end subroutine initialize_potential


    subroutine deallocate_potential(pt)
      implicit none
      TYPE(potential) :: pt
      if(associated(pt%vpot)) deallocate(pt%vpot)
      nullify(pt%vpot)
      if(associated(pt%vpotq)) deallocate(pt%vpotq)
      nullify(pt%vpotq)
      if(associated(pt%wpotq)) deallocate(pt%wpotq)
      nullify(pt%wpotq)
      if(associated(pt%ijk)) deallocate(pt%ijk)
      nullify(pt%ijk)
      return
    end subroutine deallocate_potential


    subroutine read_bands(simple_in,bd)

!this subroutine read in the bands structure from disk
!and distribute it among processors

      USE input_simple_exc, ONLY : input_options 
      USE mp,                   ONLY : mp_bcast
      USE mp_world,             ONLY : world_comm
      USE io_files,  ONLY : tmp_dir
      USE io_global, ONLY : ionode_id, ionode, stdout
      USE mp_world,  ONLY : mpime, nproc
      USE constants, ONLY : rytoev

      implicit none 
      INTEGER, EXTERNAL :: find_free_unit
      
      
      TYPE(input_options) :: simple_in
      TYPE(bands) :: bd

      INTEGER :: iun,ik,i
      INTEGER :: l_blk
      REAL(kind=DP) :: xk(3)
      REAL(kind=DP), allocatable :: et(:)
      COMPLEX(kind=DP), allocatable :: omat(:,:)

      if(ionode) then
         iun = find_free_unit()
         open( unit=iun, file=trim(tmp_dir)//trim(simple_in%prefix)//'.wfc_basis', status='old',form='unformatted')
         read(iun) bd%nk
         read(iun) bd%numv
         read(iun) bd%numc
         read(iun) bd%ntot_e
      end if

      call mp_bcast(bd%nk, ionode_id, world_comm)
      call mp_bcast(bd%numv, ionode_id, world_comm)
      call mp_bcast(bd%numc, ionode_id, world_comm)
      call mp_bcast(bd%ntot_e, ionode_id, world_comm)

      write(stdout,*) 'NUMBER OF K POINTS : ', bd%nk
      write(stdout,*) 'NUMBER OF VALENCE STATES : ', bd%numv
      write(stdout,*) 'NUMBER OF CONDUCTION STATES : ', bd%numc
      write(stdout,*) 'NUMBER OF GLOBAL STATES : ', bd%ntot_e

      bd%num=bd%numv+bd%numc

      l_blk=bd%nk/nproc
      if(l_blk*nproc<bd%nk) l_blk=l_blk+1

      if(l_blk*mpime+1 <= bd%nk) then
         bd%ik_first=l_blk*mpime+1 
         bd%ik_last=bd%ik_first+l_blk-1
         if(bd%ik_last>bd%nk) bd%ik_last=bd%nk
         bd%nk_loc=bd%ik_last-bd%ik_first+1

      else
         bd%nk_loc=0
         bd%ik_first=0
         bd%ik_last=-1
      endif
!      write(stdout,*) 'DEBUG nk_loc :', bd%nk_loc
!nk_loc ik_first ik_last
      if(bd%nk_loc>0) then
         allocate(bd%k(3,bd%nk_loc))
         allocate(bd%omat(bd%ntot_e,bd%num,bd%nk_loc))
         allocate(bd%en_v(bd%numv,bd%nk_loc))
         allocate(bd%en_c(bd%numc,bd%nk_loc))
      else
         nullify(bd%k)
         nullify(bd%omat)
         nullify(bd%en_v)
         nullify(bd%en_c)
      endif

      allocate(et(bd%num))
      allocate(omat(bd%ntot_e,bd%num))

      do ik=1,bd%nk
         if(ionode) then
            read(iun) xk(1:3)
            read(iun) et(1:bd%num)
            do i=1,bd%num
               read(iun) omat(1:bd%ntot_e,i)
            enddo
         endif
         call mp_bcast(xk, ionode_id, world_comm )
         call mp_bcast(et, ionode_id, world_comm )
         call mp_bcast(omat, ionode_id, world_comm )
         if(ik>=bd%ik_first .and. ik <= bd%ik_last) then
            bd%k(1:3,ik-bd%ik_first+1)=xk(1:3)
            bd%en_v(1:bd%numv,ik-bd%ik_first+1)=et(1:bd%numv)
            bd%en_c(1:bd%numc,ik-bd%ik_first+1)=et(bd%numv+1:bd%num)
            bd%omat(1:bd%ntot_e,1:bd%num,ik-bd%ik_first+1)=omat(1:bd%ntot_e,1:bd%num)
         endif
      enddo

      if(ionode) close(iun)

      deallocate(et,omat)

      bd%scissor=simple_in%scissor/rytoev

      return
    end subroutine read_bands

    subroutine read_product(simple_in,pd)
!read product terms from disk

      USE input_simple_exc, ONLY : input_options
      USE mp,                   ONLY : mp_bcast
      USE mp_world,             ONLY : world_comm
      USE io_files,  ONLY : tmp_dir
      USE io_global, ONLY : ionode_id, ionode, stdout
      USE mp_world,  ONLY : mpime, nproc

      implicit none

      TYPE(input_options) :: simple_in
      TYPE(product) :: pd

      INTEGER, EXTERNAL :: find_free_unit
      INTEGER :: iun
      INTEGER :: ii,jj

      if(ionode) then
         iun = find_free_unit()
         open( unit=iun, file=trim(tmp_dir)//trim(simple_in%prefix)//'.product_basis', status='old',form='unformatted')
         read(iun) pd%nprod_e
         read(iun) pd%ntot_e
      end if

      call mp_bcast(pd%nprod_e, ionode_id, world_comm)
      call mp_bcast(pd%ntot_e, ionode_id, world_comm)
    
      write(stdout,*) 'NUMBER OF PRODUCTS : ', pd%nprod_e
      write(stdout,*) 'NUMBER OF GLOBAL STATES : ', pd%ntot_e
      allocate(pd%fij(pd%nprod_e,pd%ntot_e,pd%ntot_e))
!note the order important for index corresponding to complex conjugate
      if(ionode) then
         do ii=1,pd%ntot_e
            do jj=1,pd%ntot_e
               read(iun) pd%fij(1:pd%nprod_e,ii,jj)
            enddo
         enddo
         close(iun)
      endif
      call mp_bcast(pd%fij, ionode_id, world_comm)

      return
    end subroutine read_product

    subroutine read_potential(simple_in,pt)
      
      USE input_simple_exc, ONLY : input_options
      USE mp,                   ONLY : mp_bcast
      USE mp_world,             ONLY : world_comm
      USE io_files,  ONLY : tmp_dir
      USE io_global, ONLY : ionode_id, ionode, stdout
      USE mp_world,  ONLY : mpime, nproc

      implicit none
      INTEGER, EXTERNAL :: find_free_unit


      TYPE(input_options) :: simple_in
      TYPE(potential) :: pt

      INTEGER :: iun,ii,ik,jk,kk
      INTEGER :: nktot
     

      if(ionode) then
         iun = find_free_unit()
         open( unit=iun, file=trim(tmp_dir)//trim(simple_in%prefix)//'.v_mat0', status='old',form='unformatted')
         read(iun) pt%nprod_e
      end if

      call mp_bcast(pt%nprod_e, ionode_id, world_comm)
      allocate(pt%vpot(pt%nprod_e,pt%nprod_e))
      if(ionode) then
         do ii=1,pt%nprod_e
            read(iun) pt%vpot(1:pt%nprod_e,ii)
         enddo
         
      endif
      call mp_bcast(pt%vpot, ionode_id, world_comm)
      if(ionode) read(iun) pt%nkp(1:3)
      call mp_bcast(pt%nkp, ionode_id, world_comm)
      pt%nk=pt%nkp(1)*pt%nkp(2)*pt%nkp(3)
      allocate(pt%ijk(3,pt%nk,pt%nk),pt%vpotq(pt%nprod_e,pt%nprod_e,pt%nkp(1),pt%nkp(2),pt%nkp(3)))
      allocate(pt%wpotq(pt%nprod_e,pt%nprod_e,pt%nkp(1),pt%nkp(2),pt%nkp(3)))
      if(ionode) then
         do ik=1,pt%nk!on k'
            do jk=1,pt%nk! on k
               read(iun) pt%ijk(1:3,jk,ik)
            enddo
         enddo
      endif
      call mp_bcast(pt%ijk, ionode_id, world_comm)
      do ik=1,pt%nkp(1)
         do jk=1,pt%nkp(2)
            do kk=1,pt%nkp(3)
               if(ionode) then
                  do ii=1,pt%nprod_e
                     read(iun) pt%vpotq(1:pt%nprod_e,ii,ik,jk,kk)
                  enddo
               endif
               call mp_bcast(pt%vpotq(:,:,ik,jk,kk), ionode_id, world_comm)
            enddo
         enddo
      enddo

       do ik=1,pt%nkp(1)
         do jk=1,pt%nkp(2)
            do kk=1,pt%nkp(3)
               if(ionode) then
                  do ii=1,pt%nprod_e
                     read(iun) pt%wpotq(1:pt%nprod_e,ii,ik,jk,kk)
                  enddo
               endif
               call mp_bcast(pt%wpotq(:,:,ik,jk,kk), ionode_id, world_comm)
            enddo
         enddo
      enddo


      if(ionode) close(iun)
      return
    end subroutine read_potential
    subroutine setup_exc(bd,a)
!from bands object sets up exciton a
       implicit none
       
       TYPE(bands), INTENT(in) :: bd
       TYPE(exc), INTENT(out) :: a
       a%numv=bd%numv
       a%numc=bd%numc
       a%num=bd%num
       a%nk=bd%nk
       a%nk_loc=bd%nk_loc
       a%ik_first=bd%ik_first
       a%ik_last=bd%ik_last

       if(a%nk_loc>0) then
          allocate(a%avc(a%numv,a%numc,a%nk_loc))
       else
          nullify(a%avc)
       endif
       return
    end subroutine setup_exc

    FUNCTION sum_exc(a,b)
      USE io_global, ONLY : stdout
! a+b
      implicit none
      TYPE(exc), INTENT(in) :: a,b
      TYPE(exc) :: sum_exc
!check for consistency 
      if(a%numv/=b%numv .or. a%numc/=b%numc.or. a%num/=b%num .or. a%nk/=b%nk &
           &.or.a%nk_loc/=b%nk_loc .or. a%ik_first/=b%ik_first .or. a%ik_last/=b%ik_last) then
         write(stdout,*) 'Problem with sum_exc: inconsistency'
         stop
      endif

      sum_exc%numv=b%numv
      sum_exc%numc=b%numc
      sum_exc%num=b%num
      sum_exc%nk=b%nk
      sum_exc%nk_loc=b%nk_loc
      sum_exc%ik_first=b%ik_first
      sum_exc%ik_last=b%ik_last
      if(associated(sum_exc%avc)) deallocate(sum_exc%avc)
      if(sum_exc%nk_loc>0) then
         allocate(sum_exc%avc(sum_exc%numv,sum_exc%numc,sum_exc%nk_loc))
      else
         nullify(sum_exc%avc)
      endif
      if(a%nk_loc>0) then
         sum_exc%avc(1:a%numv,1:a%numc,1:a%nk_loc)=a%avc(1:a%numv,1:a%numc,1:a%nk_loc)+&
              &b%avc(1:a%numv,1:a%numc,1:a%nk_loc)
      endif
      
      return
    END FUNCTION sum_exc
    SUBROUTINE  sum_exc_sub(sum_exc, a,b)
     USE io_global, ONLY : stdout
! a+b                                                                                                                                                         
      implicit none
      TYPE(exc), INTENT(in) :: a,b
      TYPE(exc) :: sum_exc
!check for consistency                                                                                                                                        
      if(a%numv/=b%numv .or. a%numc/=b%numc.or. a%num/=b%num .or. a%nk/=b%nk &
           &.or.a%nk_loc/=b%nk_loc .or. a%ik_first/=b%ik_first .or. a%ik_last/=b%ik_last) then
         write(stdout,*) 'Problem with sum_exc: inconsistency'
         stop
      endif

      sum_exc%numv=b%numv
      sum_exc%numc=b%numc
      sum_exc%num=b%num
      sum_exc%nk=b%nk
      sum_exc%nk_loc=b%nk_loc
      sum_exc%ik_first=b%ik_first
      sum_exc%ik_last=b%ik_last
      if(associated(sum_exc%avc)) deallocate(sum_exc%avc)
      if(sum_exc%nk_loc>0) then
         allocate(sum_exc%avc(sum_exc%numv,sum_exc%numc,sum_exc%nk_loc))
      else
         nullify(sum_exc%avc)
      endif
      if(a%nk_loc>0) then
         sum_exc%avc(1:a%numv,1:a%numc,1:a%nk_loc)=a%avc(1:a%numv,1:a%numc,1:a%nk_loc)+&
              &b%avc(1:a%numv,1:a%numc,1:a%nk_loc)
      endif

      return
    END SUBROUTINE sum_exc_sub



    FUNCTION prod_exc(a,b)
      USE io_global, ONLY : stdout
      USE mp,                   ONLY : mp_sum
      USE mp_world,             ONLY : world_comm

! scalar produc                                                                                                                                                                                          
      implicit none
      TYPE(exc), INTENT(in) :: a,b
      COMPLEX(kind=DP) :: prod_exc
      COMPLEX(kind=DP), EXTERNAL :: zdotc
!check for consistency                                                                                                                                                                         
      if(a%numv/=b%numv .or. a%numc/=b%numc.or. a%num/=b%num .or. a%nk/=b%nk &
           &.or.a%nk_loc/=b%nk_loc .or. a%ik_first/=b%ik_first .or. a%ik_last/=b%ik_last) then
         write(*,*) 'Problem with prod_exc: inconsistency'
         stop
      endif

      if(a%nk_loc>0) then
         prod_exc= zdotc(a%numv*a%numc*a%nk_loc,a%avc,1,b%avc,1)
      else
         prod_exc=(0.d0,0.d0)
      endif
      call mp_sum(prod_exc,world_comm)


      return
    END FUNCTION prod_exc
    
    FUNCTION prod_c_exc(z,a)
!complex * exciton
      USE io_global, ONLY : stdout
      USE mp,                   ONLY : mp_sum
      USE mp_world,             ONLY : world_comm

      implicit none
      TYPE(exc), INTENT(in) :: a
      COMPLEX(kind=DP), INTENT(in) :: z
      TYPE(exc) :: prod_c_exc


      prod_c_exc%numv=a%numv
      prod_c_exc%numc=a%numc
      prod_c_exc%num=a%num
      prod_c_exc%nk=a%nk
      prod_c_exc%nk_loc=a%nk_loc
      prod_c_exc%ik_first=a%ik_first
      prod_c_exc%ik_last=a%ik_last




      !if(associated(prod_c_exc%avc)) deallocate(prod_c_exc%avc)
      if(prod_c_exc%nk_loc>0) then
         allocate(prod_c_exc%avc(prod_c_exc%numv,prod_c_exc%numc,prod_c_exc%nk_loc))
      else
         nullify(prod_c_exc%avc)
      endif
      if(prod_c_exc%nk_loc>0) then
         prod_c_exc%avc(1:a%numv,1:a%numc,1:a%nk_loc)=z*a%avc(1:a%numv,1:a%numc,1:a%nk_loc)

      endif


      return
    END FUNCTION prod_c_exc

    SUBROUTINE assign_exc(a,b)
!x=a
      implicit none

      TYPE(exc), INTENT(out) :: a
      TYPE(exc),INTENT(in) ::b


      a%numv=b%numv
      a%numc=b%numc
      a%num=b%num
      a%nk=b%nk
      a%nk_loc=b%nk_loc
      a%ik_first=b%ik_first
      a%ik_last=b%ik_last




      if(associated(a%avc)) deallocate(a%avc)
      if(a%nk_loc>0) then
         allocate(a%avc(a%numv,a%numc,a%nk_loc))
      else
         nullify(a%avc)
      endif
      if(a%nk_loc>0) then
         a%avc(1:a%numv,1:a%numc,1:a%nk_loc)=b%avc(1:a%numv,1:a%numc,1:a%nk_loc)

      endif

      return
    END SUBROUTINE assign_exc

    SUBROUTINE randomize_exc(a)
!this subroutine set an exc object randomly

      USE random_numbers, ONLY : randy

      implicit none
      TYPE(exc), INTENT(inout) :: a
      INTEGER :: i,j,k

      if(a%nk_loc > 0) then
         if(a%nk_loc >  0 ) then
            do i=1,a%numv
               do j=1,a%numc
                  do k=1,a%nk_loc
                     a%avc(i,j,k)=dcmplx(randy(),randy())
                  enddo
               enddo
            enddo
         endif
      endif
      return
    END SUBROUTINE randomize_exc

    SUBROUTINE normalize_exc(a)
!normalize the avc vector
      implicit none
      TYPE(exc), INTENT(inout) :: a

      COMPLEX(kind=DP) :: csca
      REAL(kind=DP) :: sca

      csca=a*a
      sca=dble(csca)
      csca=cmplx(1.d0/sqrt(sca),0.d0)

      if(a%nk_loc>0)    a%avc=csca*a%avc
      return

    END SUBROUTINE normalize_exc

    FUNCTION h_diag(bd,a)
!this function applies the diagonal part of the excitonic Hamiltonian
!TO DO SCISSOR OR SIMILA
      USE io_global, ONLY : stdout,ionode

      implicit none
      TYPE(exc) :: h_diag
      TYPE(exc), INTENT(in) :: a
      TYPE(bands), INTENT(in) :: bd
      INTEGER :: iv,ic,ik

      h_diag%numv=a%numv
      h_diag%numc=a%numc
      h_diag%num=a%num
      h_diag%nk=a%nk
      h_diag%nk_loc=a%nk_loc
      h_diag%ik_first=a%ik_first
      h_diag%ik_last=a%ik_last



      
      if(h_diag%nk_loc>0) then
         !if(associated(h_diag%avc)) deallocate(h_diag%avc)
         allocate(h_diag%avc(h_diag%numv,h_diag%numc,h_diag%nk_loc))
      else
         nullify(h_diag%avc)
      endif
      if(ionode) then
 !        write(stdout,*) 'DEBUG:',h_diag%nk_loc,h_diag%numc,h_diag%numv,bd%en_c(1,1)
 !        write(stdout,*) 'DEBUG:',bd%nk_loc,bd%numc,bd%numv,bd%en_v(1,1)
      endif
      
      if(h_diag%nk_loc>0) then
         do ik=1,h_diag%nk_loc
            do ic=1,h_diag%numc
                do iv=1,h_diag%numv
                   h_diag%avc(iv,ic,ik)=(bd%en_c(ic,ik)-bd%en_v(iv,ik)+bd%scissor)*a%avc(iv,ic,ik)  
                enddo
             enddo
          enddo

      endif



      return
    END FUNCTION h_diag
END MODULE simple_objects
