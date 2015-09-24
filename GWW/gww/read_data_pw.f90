!
! Copyright (C) 2001-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!



!these subroutines read in the data from PW calculations



   SUBROUTINE read_data_pw_u(wu,prefix)
!this subroutine reads in the energies and the inversewannier transformation matrix

    USE kinds,                ONLY : DP
    USE io_global,            ONLY : stdout, ionode, ionode_id
    USE basic_structures
    USE mp,                   ONLY : mp_bcast, mp_barrier
    USE mp_world,             ONLY : world_comm
    USE io_files,  ONLY : tmp_dir

    implicit none
    INTEGER, EXTERNAL :: find_free_unit
    TYPE(wannier_u) :: wu!structure to be read and initialized
    CHARACTER(LEN=256) ::  prefix!to designate the PW files

    INTEGER :: iunu
    INTEGER :: iw,is

    if(ionode) then
       iunu = find_free_unit()
       open( unit=iunu, file=trim(tmp_dir)//trim(prefix)//'.wannier', status='old',form='unformatted')

!read in basis length
       read(iunu) wu%nspin
       read(iunu) wu%nums
    endif

    call mp_bcast(wu%nspin, ionode_id, world_comm)
    call mp_bcast(wu%nums, ionode_id, world_comm)


!allocate arrays    
    allocate(wu%ene(wu%nums,wu%nspin))
    allocate(wu%ene_xc(wu%nums,wu%nspin))
    allocate(wu%ene_lda_h(wu%nums,wu%nspin))
    allocate(wu%umat(wu%nums,wu%nums,wu%nspin))



    do is=1,wu%nspin
       if(ionode)   read(iunu) wu%nums_occ(is)
       !write(stdout,*) 'DEBUG:', wu%nspin,wu%nums,wu%nums_occ(is)
       !FLUSH(stdout)

       call mp_bcast(wu%nums_occ(is), ionode_id, world_comm)
      
       if(ionode) then
!read in energies
          read(iunu) wu%ene(1:wu%nums,is)
!read in DFT exchange and correlation energies

          read(iunu) wu%ene_xc(1:wu%nums,is)
          read(iunu) wu%ene_lda_h(1:wu%nums,is)
!read in transformation matrix
          do iw=1,wu%nums
             read(iunu) wu%umat(1:wu%nums,iw,is)
          enddo
       endif

       call mp_bcast(wu%ene(:,is), ionode_id, world_comm)
       call mp_bcast(wu%ene_xc(:,is), ionode_id, world_comm)
       call mp_bcast(wu%ene_lda_h(:,is), ionode_id, world_comm)
       do iw=1,wu%nums
          call mp_barrier( world_comm )
          call mp_bcast(wu%umat(:,iw,is), ionode_id, world_comm)
       enddo
    enddo
    if(ionode) close(iunu)
    return

  END

  SUBROUTINE  read_data_pw_v(vp,prefix,debug,ort,l_zero)
!read from file and initialize coulomb potential on the basis of products of wanniers

    USE kinds,                ONLY : DP
    USE io_global,            ONLY : stdout, ionode, ionode_id
    USE constants,            ONLY : eps8
    USE basic_structures
    USE mp,                   ONLY : mp_bcast,mp_barrier
    USE mp_world,             ONLY : world_comm
   USE io_files,  ONLY :  tmp_dir

    implicit none
    INTEGER, EXTERNAL :: find_free_unit
    TYPE(v_pot) :: vp!structure to be read and initialized
    CHARACTER(LEN=256) ::  prefix!to designate the PW files
    LOGICAL :: debug!if true check for simmetry
    INTEGER :: ort!if ort==0, open non orthogonal file, if ort == 1 open orthogonal file
                  !if ort==2 open non orthogonal symmetric file
    LOGICAL :: l_zero!if true open file with head put to zero

    INTEGER :: iunv
    INTEGER :: iw,jw

    if(ionode) then
       iunv=find_free_unit()
       if(ort==1) then
          open( unit=iunv, file=trim(tmp_dir)//trim(prefix)//'.vpot', status='old',form='unformatted')
       else if (ort==0) then
          if(.not.l_zero) then
             open( unit=iunv, file=trim(tmp_dir)//trim(prefix)//'.vpot_no', status='old',form='unformatted')
          else
             open( unit=iunv, file=trim(tmp_dir)//trim(prefix)//'.vpot_no_zero', status='old',form='unformatted')
          endif
       else if (ort==2) then
          if(.not.l_zero) then
             open( unit=iunv, file=trim(tmp_dir)//trim(prefix)//'.vpot_no_sym', status='old',form='unformatted')
          else
             open( unit=iunv, file=trim(tmp_dir)//trim(prefix)//'.vpot_no_sym_zero', status='old',form='unformatted')
          endif
       endif
!read in basis length
       read(iunv) vp%numpw

    endif

    call mp_bcast(vp%numpw, ionode_id, world_comm)

!allocate array
    allocate(vp%vmat(vp%numpw,vp%numpw))


!read in potential matrix
    if(ionode) then
       do iw=1,vp%numpw
          read(iunv) vp%vmat(1:vp%numpw,iw)
       enddo
    endif

    do iw=1,vp%numpw
       call mp_barrier( world_comm )
       call mp_bcast(vp%vmat(:,iw), ionode_id, world_comm)
    enddo
!check
    if(debug) then
      do iw=1,vp%numpw
        do jw=1,iw
           if(abs(vp%vmat(iw,jw)-vp%vmat(jw,iw)) >= eps8) then
             write(stdout,*) 'Proble vmat not simmetric:',iw,jw,vp%vmat(iw,jw)-vp%vmat(jw,iw)
           endif
        enddo
      enddo
    endif

    if(ionode) close(iunv)

    return
 END SUBROUTINE
 
 SUBROUTINE read_data_pw_q(qm,prefix,l_v_products)
!this subroutine reads in and allocate the arrays for the 
!description of overlaps of (orthonormalized) products of wanniers
!with products of wannier


    USE kinds,                ONLY : DP
    USE io_global,            ONLY : stdout, ionode, ionode_id
    USE basic_structures
    USE mp,                   ONLY : mp_bcast, mp_barrier
    USE mp_world,             ONLY : world_comm
   USE io_files,  ONLY :  tmp_dir

    implicit none
    INTEGER, EXTERNAL :: find_free_unit
    TYPE(q_mat) :: qm!structure to be read and initialized
    CHARACTER(LEN=256) ::  prefix!to designate the PW files
    LOGICAL :: l_v_products!if true read the wp_v file for the products \tilde{w}_i(r)\tilde{w}_j(r)v(r,r')w^P_red_k(r')

    INTEGER :: iunq
    INTEGER :: iw


    if(ionode) then
       iunq=find_free_unit()
       if(.not.l_v_products) then
          open( unit=iunq, file=trim(tmp_dir)//trim(prefix)//'.wp', status='old',form='unformatted')
       else
          open( unit=iunq, file=trim(tmp_dir)//trim(prefix)//'.wp_v', status='old',form='unformatted')
       endif

!read in basis length
       read(iunq) qm%numpw
    endif

    call mp_bcast(qm%numpw, ionode_id, world_comm)

! allocate array of descriptors
    allocate (qm%wp(qm%numpw))


    do iw=1,qm%numpw

       if(ionode)    read(iunq) qm%wp(iw)%numij
       call mp_bcast(qm%wp(iw)%numij, ionode_id, world_comm)
       
!for each descriptor allocates arrays
       allocate(qm%wp(iw)%ij(2,qm%wp(iw)%numij))
       allocate(qm%wp(iw)%o(qm%wp(iw)%numij))

!read data
       if(ionode) then
          read(iunq) qm%wp(iw)%ij(1,1:qm%wp(iw)%numij)
          read(iunq) qm%wp(iw)%ij(2,1:qm%wp(iw)%numij)
          read(iunq) qm%wp(iw)%o(1:qm%wp(iw)%numij)
       end if
       call mp_bcast(qm%wp(iw)%ij(:,:), ionode_id, world_comm)
       call mp_bcast(qm%wp(iw)%o(:), ionode_id, world_comm)
    enddo


    qm%is_parallel=.false.
    qm%numpw_para=qm%numpw
    qm%first_para=1

    if(ionode) close(iunq)

    return

 END SUBROUTINE

 SUBROUTINE read_data_pw_ortho_polaw(op,prefix)
!this subroutine reads in and allocate the arrays for the
!description of orthonormalization matrix of wanniers products


    USE kinds,                ONLY : DP
    USE io_global,            ONLY : stdout, ionode, ionode_id
    USE basic_structures,     ONLY : ortho_polaw, free_memory
    USE mp,                   ONLY : mp_bcast, mp_barrier
    USE mp_world,             ONLY : world_comm
   USE io_files,  ONLY :  tmp_dir

    implicit none
    INTEGER, EXTERNAL :: find_free_unit
    TYPE(ortho_polaw) :: op!structure to be read and initialized
    CHARACTER(LEN=256) ::  prefix!to designate the PW files

    INTEGER :: iunq
    INTEGER :: iw

!    call free_memory(op)

    if(ionode) then
       iunq=find_free_unit()
       open( unit=iunq, file=trim(tmp_dir)//trim(prefix)//'.orthonorm', status='old',form='unformatted')

!read in basis length
       read(iunq) op%numpw
    endif
    call mp_bcast(op%numpw, ionode_id, world_comm)
    allocate(op%on_mat(op%numpw,op%numpw))

    if(ionode) then
       do iw=1,op%numpw
          read(iunq) op%on_mat(1:op%numpw,iw)
       enddo
    end if
    do iw=1,op%numpw
       call mp_barrier( world_comm )
       call mp_bcast( op%on_mat(:,iw), ionode_id, world_comm)
    enddo

    op%inverse=.false.

    if(ionode) close(iunq)

 END subroutine

 SUBROUTINE read_data_pw_wp_psi(wp,prefix)
!this subroutine reads in and allocate the arrays for the
!description of products of valence^2 times two wannier products 

    USE kinds,                ONLY : DP
    USE io_global,            ONLY : stdout, ionode, ionode_id
    USE basic_structures,     ONLY : wp_psi, free_memory
    USE mp,                   ONLY : mp_bcast
    USE mp_world,             ONLY : world_comm
   USE io_files,  ONLY : tmp_dir

    implicit none
    INTEGER, EXTERNAL :: find_free_unit
    TYPE(wp_psi) :: wp!structure to be read and initialized
    CHARACTER(LEN=256) ::  prefix!to designate the PW files

    INTEGER :: iunq
    INTEGER :: iw,hw, jw


!    call free_memory(wp)

    if(ionode) then
       iunq=find_free_unit()
       open( unit=iunq, file=trim(tmp_dir)//trim(prefix)//'.wpwp_psi', status='old',form='unformatted')

!read in basis length
       read(iunq) wp%numpw
       read(iunq) wp%nums_psi     
    endif
    call mp_bcast(wp%numpw, ionode_id, world_comm)
    call mp_bcast(wp%nums_psi, ionode_id, world_comm)
    allocate(wp%wwp(wp%numpw,wp%numpw,wp%nums_psi))

    do hw=1,wp%nums_psi
       if(ionode) then
          do iw=1,wp%numpw
             read(iunq) wp%wwp(iw,1:iw,hw)
          enddo
          do iw=1,wp%numpw
             do jw=iw, wp%numpw
                 wp%wwp(iw,jw,hw)=wp%wwp(jw,iw,hw)
              enddo
           enddo
        endif
        call mp_bcast( wp%wwp(:,:,hw), ionode_id, world_comm)
     enddo

    if(ionode) close(iunq)

 END subroutine


   SUBROUTINE read_data_pw_u_prim(wu,prefix)
!this subroutine reads in the inverse wannier transformation matrix

    USE kinds,                ONLY : DP
    USE io_global,            ONLY : stdout, ionode, ionode_id
    USE basic_structures
    USE mp,                   ONLY : mp_bcast, mp_barrier
    USE mp_world,             ONLY : world_comm
   USE io_files,  ONLY :  tmp_dir

    implicit none
    INTEGER, EXTERNAL :: find_free_unit
    TYPE(wannier_u_prim) :: wu!structure to be read and initialized
    CHARACTER(LEN=256) ::  prefix!to designate the PW files

    INTEGER :: iunu
    INTEGER :: iw

    if(ionode) then
       iunu = find_free_unit()
       open( unit=iunu, file=trim(tmp_dir)//trim(prefix)//'.wannier_prim', status='old',form='unformatted')

!read in basis length
       read(iunu) wu%nums_prim
       read(iunu) wu%nums_occ
       read(iunu) wu%nums
       write(*,*) 'read_data_pw_u_prim',wu%nums_prim,wu%nums_occ,wu%nums
    endif

    call mp_bcast(wu%nums_prim, ionode_id, world_comm)
    call mp_bcast(wu%nums_occ, ionode_id, world_comm)
    call mp_bcast(wu%nums, ionode_id, world_comm)

!allocate arrays
    allocate(wu%umat(wu%nums_prim,wu%nums_prim))

    if(ionode) then
!read in transformation matrix
       do iw=1,wu%nums_prim
          read(iunu) wu%umat(1:wu%nums_prim,iw)
       enddo
    endif

    do iw=1,wu%nums_prim
       call mp_barrier( world_comm )
       call mp_bcast(wu%umat(:,iw), ionode_id, world_comm)
    enddo

    if(ionode) close(iunu)
    return

  END SUBROUTINE read_data_pw_u_prim


   SUBROUTINE read_data_pw_v_pot_prim(vp,prefix, l_zero)
!this subroutine reads in the coulombian potential and the overlap index 

    USE kinds,                ONLY : DP
    USE io_global,            ONLY : stdout, ionode, ionode_id
    USE basic_structures
    USE mp,                   ONLY : mp_bcast, mp_barrier
    USE mp_world,             ONLY : world_comm
   USE io_files,  ONLY :  tmp_dir

    implicit none
    INTEGER, EXTERNAL :: find_free_unit
    TYPE(v_pot_prim) :: vp!structure to be read and initialized
    CHARACTER(LEN=256) ::  prefix!to designate the PW files
    LOGICAL :: l_zero!if true opens file with head of v pu to zero

    INTEGER :: iunu
    INTEGER :: iw

    if(ionode) then
       iunu = find_free_unit()
       if(.not. l_zero) then
          open( unit=iunu, file=trim(tmp_dir)//trim(prefix)//'.uterms_prim', status='old',form='unformatted')
       else
          open( unit=iunu, file=trim(tmp_dir)//trim(prefix)//'.uterms_prim_zero', status='old',form='unformatted')
       endif

!read in basis length
       read(iunu) vp%numpw_prim
       read(iunu) vp%numpw
       write(*,*) 'read_data_pw_v_pot_prim', vp%numpw_prim,vp%numpw
    endif

    call mp_bcast(vp%numpw, ionode_id, world_comm)
    call mp_bcast(vp%numpw_prim, ionode_id, world_comm)

!allocate arrays
    allocate(vp%vmat(vp%numpw_prim,vp%numpw))
    allocate(vp%ij(2,vp%numpw_prim))

    if(ionode) then
!read in transformation matrix
       do iw=1,vp%numpw_prim
          read(iunu) vp%vmat(iw,1:vp%numpw)
       enddo
       close(iunu)
    endif

    do iw=1,vp%numpw
       call mp_barrier( world_comm )
       call mp_bcast(vp%vmat(:,iw), ionode_id, world_comm)
    enddo

    if(ionode) then
       iunu = find_free_unit()
       open( unit=iunu, file=trim(tmp_dir)//trim(prefix)//'.ij_prim', status='old',form='unformatted')
       do iw=1,vp%numpw_prim
          read(iunu) vp%ij(1,iw),vp%ij(2,iw)
       enddo
       close(iunu)
    endif


    call mp_bcast(vp%ij(:,:), ionode_id, world_comm)

    vp%is_parallel=.false.
    vp%numpw_para=vp%numpw
    vp%first_para=1


    return

  END SUBROUTINE read_data_pw_v_pot_prim

  SUBROUTINE read_data_pw_wp_psi_cutoff_index(wpi,prefix)
!this subroutine reads in and allocate the arrays for the
!indices describing of products of valence^2 times two wannier products
!when a cutoff is applied

    USE kinds,                ONLY : DP
    USE io_global,            ONLY : stdout, ionode, ionode_id
    USE basic_structures,     ONLY : wp_psi_cutoff_index, free_memory
    USE mp,                   ONLY : mp_bcast
    USE mp_world,             ONLY : world_comm
   USE io_files,  ONLY :  tmp_dir

    implicit none
    INTEGER, EXTERNAL :: find_free_unit
    TYPE(wp_psi_cutoff_index) :: wpi!structure to be read and initialized
    CHARACTER(LEN=256) ::  prefix!to designate the PW files

    INTEGER :: iuni
    INTEGER :: i

    if(ionode) then
       iuni=find_free_unit()
       open( unit=iuni, file=trim(tmp_dir)//trim(prefix)//'.wpwp_psi_index', status='old',form='unformatted')
!read in basis length
       read(iuni) wpi%numpw
       read(iuni) wpi%nums_psi
       read(iuni) wpi%numpwpw
    endif

    call mp_bcast(wpi%numpw, ionode_id, world_comm)
    call mp_bcast(wpi%nums_psi, ionode_id, world_comm)
    call mp_bcast(wpi%numpwpw, ionode_id, world_comm)
    
    allocate(wpi%index(2,wpi%numpwpw))

    if(ionode) then
       do i=1,wpi%numpwpw
          read(iuni) wpi%index(1,i),wpi%index(2,i)
       enddo
       close(iuni)
    endif

    call mp_bcast(wpi%index, ionode_id, world_comm)
    
    return
  END SUBROUTINE read_data_pw_wp_psi_cutoff_index

  SUBROUTINE read_data_pw_wp_psi_cutoff_data(wpi,wp,prefix)
!this subroutine reads in and allocate the arrays for the
!products of valence^2 times two wannier products when a cutoff is applied

    USE kinds,                ONLY : DP
    USE io_global,            ONLY : stdout, ionode, ionode_id
    USE basic_structures,     ONLY : wp_psi_cutoff_index, wp_psi_cutoff_data,free_memory
    USE mp,                   ONLY : mp_bcast
    USE mp_world,             ONLY : world_comm
   USE io_files,  ONLY :  tmp_dir

    implicit none
    INTEGER, EXTERNAL :: find_free_unit
    TYPE(wp_psi_cutoff_index), INTENT(in) :: wpi!indices
    TYPE(wp_psi_cutoff_data), INTENT(inout) :: wp!data to be read
    CHARACTER(LEN=256), INTENT(in) ::  prefix!to designate the PW files

    INTEGER :: iund
    INTEGER :: i, pos,state
    REAL(kind=DP) :: w

    
    wp%numpw=wpi%numpw
    wp%nums_psi=wpi%nums_psi
    wp%numpwpw=wpi%numpwpw
    
    allocate(wp%wwp(wp%numpwpw,wp%nums_psi))

    if(ionode) then
       iund=find_free_unit()
       open( unit=iund, file=trim(tmp_dir)//trim(prefix)//'.wpwp_psi', status='old',form='unformatted')
!read in basis length
       do i=1,wp%nums_psi*wp%numpwpw
          read(iund) pos,state,w
          wp%wwp(pos,state)=w
       enddo
       close(iund)
    endif

    do i=1,wp%nums_psi
       call mp_bcast(wp%wwp(:,i), ionode_id, world_comm)
    enddo
    
    return
  END SUBROUTINE read_data_pw_wp_psi_cutoff_data
    


   SUBROUTINE read_data_pw_exchange(ene_x,max_i,prefix,nspin)
!this subroutine reads in the exchange energies

    USE kinds,                ONLY : DP
    USE io_global,            ONLY : stdout, ionode, ionode_id
    USE basic_structures
    USE mp,                   ONLY : mp_bcast
    USE mp_world,             ONLY : world_comm
   USE io_files,  ONLY :  tmp_dir

    implicit none
    INTEGER, EXTERNAL :: find_free_unit
    REAL(kind=DP) :: ene_x(max_i,nspin)
    INTEGER :: max_i
    CHARACTER(LEN=256) ::  prefix!to designate the PW files
    INTEGER, INTENT(in) :: nspin!spin multiplicity

    INTEGER :: iunu
    INTEGER :: ndata,is
   
    REAL(kind=DP), ALLOCATABLE :: buf(:)

    if(ionode) then
       iunu = find_free_unit()
       open( unit=iunu, file=trim(tmp_dir)//trim(prefix)//'.exchange', status='old',form='unformatted')
       read(iunu) ndata
       allocate(buf(ndata))
       do is=1,nspin
          read(iunu) buf(1:ndata)
          ene_x(1:max_i,is)=buf(1:max_i)
       enddo
       close(iunu)
       deallocate(buf)
    endif
    call mp_bcast(ene_x, ionode_id,world_comm)


    return
  end SUBROUTINE read_data_pw_exchange


  SUBROUTINE read_data_pw_exchange_off(ene_x_off,max_i,prefix,nspin)
!this subroutine reads in the whole fock matrix                                                                                                                        

    USE kinds,                ONLY : DP
    USE io_global,            ONLY : stdout, ionode, ionode_id
    USE basic_structures
    USE mp,                   ONLY : mp_bcast
    USE mp_world,             ONLY : world_comm
   USE io_files,  ONLY :  tmp_dir

    implicit none
    INTEGER, EXTERNAL :: find_free_unit
    REAL(kind=DP) :: ene_x_off(max_i,max_i,nspin)
    INTEGER :: max_i
    CHARACTER(LEN=256) ::  prefix!to designate the PW files                                                                                                            
    INTEGER, INTENT(in) :: nspin!spin multiplicity                                                                                                                     

    INTEGER :: iunu
    INTEGER :: ndata,is,ii

    REAL(kind=DP), ALLOCATABLE :: buf(:)

    if(ionode) then
       iunu = find_free_unit()
       open( unit=iunu, file=trim(tmp_dir)//trim(prefix)//'.exchange_off', status='old',form='unformatted')
       read(iunu) ndata
       allocate(buf(ndata))
       do is=1,nspin
          do ii=1,ndata
             read(iunu) buf(1:ndata)
             if(ii<=max_i) ene_x_off(1:max_i,ii,is)=buf(1:max_i)
          enddo
       enddo
       close(iunu)
       deallocate(buf)
    endif
    call mp_bcast(ene_x_off, ionode_id,world_comm)


    return
  end SUBROUTINE read_data_pw_exchange_off




  SUBROUTINE read_data_pw_head_epsilon(he, prefix, l_wing_epsilon, l_gzero)
!this subroutine reads the data


    USE kinds,                ONLY : DP
    USE io_global,            ONLY : stdout, ionode, ionode_id
    USE basic_structures,     ONLY : head_epsilon
    USE mp,                   ONLY : mp_bcast, mp_barrier
    USE mp_world,             ONLY : world_comm
   USE io_files,  ONLY :  tmp_dir

    implicit none
    INTEGER, EXTERNAL :: find_free_unit
    TYPE(head_epsilon) :: he!the head of epsilon to be read
    CHARACTER(LEN=256) ::  prefix!to designate the PW files
    LOGICAL :: l_wing_epsilon!if true read from file also the wing data
    LOGICAL :: l_gzero!if true reads also gzero otherwise is initialized to 0

    INTEGER :: iun,i, idumm,ipol
    REAL(kind=DP) :: rdumm

    if(ionode) then
       iun = find_free_unit()
       open( unit=iun, file=trim(tmp_dir)//'/_ph0/'//trim(prefix)//'.head', status='old',form='unformatted')
       read(iun) he%n
       read(iun) he%omega
    endif
    call mp_bcast(he%n, ionode_id,world_comm)
    call mp_bcast(he%omega, ionode_id,world_comm)
    allocate(he%freqs(he%n+1))
    allocate(he%head(he%n+1,3))
    if(ionode) then
       read(iun) he%freqs(1:he%n+1)
       do ipol=1,3
          read(iun) he%head(1:he%n+1,ipol)
       enddo
       close(iun)
    endif
    call mp_bcast(he%freqs, ionode_id,world_comm)
    call mp_bcast(he%head, ionode_id,world_comm)
    
    
    if(ionode) then
       iun = find_free_unit()
       open( unit=iun, file=trim(tmp_dir)//trim(prefix)//'.wing', status='old',form='unformatted')
       read(iun) idumm
       read(iun) rdumm
       if(idumm /= he%n) then
          write(stdout,*) 'WING: PROBLEM WITH N'
       endif
       if(rdumm /= he%omega) then
          write(stdout,*) 'WING: PROBLEM WITH OMEGA'
       endif
       read(iun) he%numpw
    endif
    call mp_bcast(he%numpw, ionode_id,world_comm)
    allocate(he%wing(he%numpw, he%n+1,3))
    allocate(he%wing_c(he%numpw, he%n+1,3))
    
!          if(idumm /= he%numpw) then
!             write(stdout,*) 'WING: PROBLEM WITH NUMPW', idumm, he%numpw
!          endif
    if(ionode) then
       do ipol=1,3
          do i=1,he%n+1
             read(iun) he%wing(1:he%numpw,i,ipol)
          enddo
       enddo
!       do i=1,he%n+1
!          read(iun) he%wing_c(1:he%numpw,i,ipol)
!       enddo
          
       close(iun)
    endif
!    do i=1,he%n+1
!       call mp_barrier
!
!       call mp_bcast(he%wing(:,i), ionode_id,world_comm)
!      call mp_bcast(he%wing_c(:,i), ionode_id,world_comm)
!    enddo
    call mp_bcast(he%wing, ionode_id,world_comm)

    if(l_gzero) then
       if(ionode) then
          iun = find_free_unit()
          open( unit=iun, file=trim(tmp_dir)//trim(prefix)//'.gzero', status='old',form='unformatted')
          read(iun) idumm
           if(idumm /= he%numpw) then
             write(stdout,*) 'WING: PROBLEM WITH NUMPW', idumm, he%numpw
          endif
       endif
       allocate(he%gzero(he%numpw))
       if(ionode) then
          do i=1,he%numpw
             read(iun) he%gzero(i)
          enddo
          close(iun)
       endif
       call mp_bcast(he%gzero,ionode_id,world_comm)
    else
       allocate(he%gzero(he%numpw))
       he%gzero(1:he%numpw)=0.d0
    endif


    return
  END SUBROUTINE read_data_pw_head_epsilon


  SUBROUTINE read_data_pw_cprim_prod(cpp, prefix, l_vc, ok_read, l_vcw_overlap, l_upper)
!this subroutine read the products cprim c v\tilde{w^P} from disk
    USE kinds,                ONLY : DP
    USE io_global,            ONLY : stdout, ionode, ionode_id
    USE basic_structures,     ONLY : cprim_prod,free_memory
    USE mp,                   ONLY : mp_bcast, mp_barrier
    USE mp_world,             ONLY : world_comm
   USE io_files,  ONLY :  tmp_dir

    implicit none

    INTEGER, EXTERNAL :: find_free_unit
    TYPE(cprim_prod) :: cpp!the structure to be read
    CHARACTER(LEN=256) ::  prefix!to designate the PW files
    LOGICAL ::  l_vc !if true reads in the vc terms
    LOGICAL, INTENT(out) :: ok_read!if true effectively read otherwise the file doesn't exist
    LOGICAL, INTENT(in) :: l_vcw_overlap!if true read the overlaps v c w
    LOGICAL, INTENT(in) :: l_upper!if true reads data for reduced upper states

    CHARACTER(4) :: nfile
    INTEGER :: iunsterms, icp, i


    call free_memory(cpp)
    if(.not.l_vcw_overlap) then
       if(ionode) then
          write(nfile,'(4i1)') &
               & cpp%cprim/1000,mod(cpp%cprim,1000)/100,mod(cpp%cprim,100)/10,mod(cpp%cprim,10)
          if(.not.l_upper) then
             if(.not. l_vc) then
                inquire(file=trim(tmp_dir)//trim(prefix)//'.cprim.'//nfile,exist=ok_read)
             else
                inquire(file=trim(tmp_dir)//trim(prefix)//'.vcprim.'//nfile,exist=ok_read)
             endif
          else
             if(.not. l_vc) then
                inquire(file=trim(tmp_dir)//trim(prefix)//'.cprim_up.'//nfile,exist=ok_read)
             else
                inquire(file=trim(tmp_dir)//trim(prefix)//'.vcprim_up.'//nfile,exist=ok_read)
             endif
          endif
       endif
       call mp_bcast(ok_read, ionode_id,world_comm)
       if(.not. ok_read) return
    endif

    if(ionode) then
       iunsterms =  find_free_unit()
       write(nfile,'(4i1)') &
            & cpp%cprim/1000,mod(cpp%cprim,1000)/100,mod(cpp%cprim,100)/10,mod(cpp%cprim,10)
       if(.not.l_upper) then
          if(l_vcw_overlap) then
             open( unit= iunsterms, file=trim(tmp_dir)//trim(prefix)//'.vcw_overlap.'//nfile, status='old',form='unformatted')
          else
             if(.not. l_vc) then
                open( unit= iunsterms, file=trim(tmp_dir)//trim(prefix)//'.cprim.'//nfile, status='old',form='unformatted')
             else
                open( unit= iunsterms, file=trim(tmp_dir)//trim(prefix)//'.vcprim.'//nfile, status='old',form='unformatted')
             endif
          endif
       else
           if(l_vcw_overlap) then
             open( unit= iunsterms, file=trim(tmp_dir)//trim(prefix)//'.vcw_up_overlap.'//nfile, status='old',form='unformatted')
          else
             if(.not. l_vc) then
                open( unit= iunsterms, file=trim(tmp_dir)//trim(prefix)//'.cprim_up.'//nfile, status='old',form='unformatted')
             else
                open( unit= iunsterms, file=trim(tmp_dir)//trim(prefix)//'.vcprim_up.'//nfile, status='old',form='unformatted')
             endif
          endif
       endif
       read(iunsterms) icp
       if(icp /= cpp%cprim) then
          write(stdout,*) 'PROBLEM WITH CPRIM_PROD'
          stop
       endif
       read(iunsterms) cpp%nums_occ
       read(iunsterms) cpp%nums!DIFFERENT MEANING FOR UPPER STATES
       read(iunsterms) cpp%numpw
    endif
    call mp_bcast(cpp%nums_occ, ionode_id,world_comm)
    call mp_bcast(cpp%nums, ionode_id,world_comm)
    call mp_bcast(cpp%numpw, ionode_id,world_comm)

    cpp%nums_cond=cpp%nums-cpp%nums_occ
    if(.not.l_vc .or. l_vcw_overlap .and. .not.l_upper) then
       allocate(cpp%cpmat(cpp%numpw,cpp%nums_cond))
    else
       allocate(cpp%cpmat(cpp%numpw,cpp%nums))
    endif
    cpp%lda=cpp%numpw
    if(.not. l_vc .or. l_vcw_overlap .and. .not.l_upper) then
       do i=1,cpp%nums_cond
          call mp_barrier( world_comm )
          if(ionode) read(iunsterms) cpp%cpmat(1:cpp%numpw,i)
          call mp_bcast(cpp%cpmat(:,i), ionode_id,world_comm)
       enddo
    else
        do i=1,cpp%nums
           call mp_barrier( world_comm )
          if(ionode) read(iunsterms) cpp%cpmat(1:cpp%numpw,i)
          call mp_bcast(cpp%cpmat(:,i), ionode_id,world_comm)
       enddo
    endif
    if(ionode) close(iunsterms)
    
    cpp%is_parallel=.false.
    cpp%numpw_para=cpp%numpw
    cpp%first_para=1

    return
  END SUBROUTINE read_data_pw_cprim_prod


  SUBROUTINE read_data_pw_dft_xc(ene_dft_xc,max_i,prefix)
!this subroutine reads in the exchange energies

    USE kinds,                ONLY : DP
    USE io_global,            ONLY : stdout, ionode, ionode_id
    USE basic_structures
    USE mp,                   ONLY : mp_bcast
    USE mp_world,             ONLY : world_comm
   USE io_files,  ONLY :  tmp_dir

    implicit none
    INTEGER, EXTERNAL :: find_free_unit
    REAL(kind=DP) :: ene_dft_xc(max_i)
    INTEGER :: max_i
    CHARACTER(LEN=256) ::  prefix!to designate the PW files

    INTEGER :: iunu, i,nn


    if(ionode) then
       iunu = find_free_unit()
       open( unit=iunu, file=trim(tmp_dir)//trim(prefix)//'.dft_xc', status='old',form='unformatted')
       read(iunu) nn
       do i=1,max_i
          read(iunu) ene_dft_xc(i)
       enddo
    endif
    call mp_bcast(ene_dft_xc(1:max_i), ionode_id, world_comm)


    return
  end SUBROUTINE read_data_pw_dft_xc


  SUBROUTINE read_data_pw_dft_xc_off(ene_dft_xc_off,max_i,prefix,ispin)
!this subroutine reads in the exchange energies                                                                                                             

    USE kinds,                ONLY : DP
    USE io_global,            ONLY : stdout, ionode, ionode_id
    USE basic_structures
    USE mp,                   ONLY : mp_bcast
    USE mp_world,             ONLY : world_comm
   USE io_files,  ONLY :  tmp_dir

    implicit none
    INTEGER, EXTERNAL :: find_free_unit
    REAL(kind=DP) :: ene_dft_xc_off(max_i,max_i)
    INTEGER :: max_i
    CHARACTER(LEN=256) ::  prefix!to designate the PW files   
    INTEGER, INTENT(in) :: ispin! spin channel

    INTEGER :: iunu, ibnd,nn


    if(ionode) then
       iunu = find_free_unit()
       if(ispin==1) then
          open(unit=iunu,file=trim(tmp_dir)//trim(prefix)//'.exc_off',status='old',form='unformatted')
       else
          open(unit=iunu,file=trim(tmp_dir)//trim(prefix)//'.exc_off2',status='old',form='unformatted')
       endif
       read(iunu) nn
       do ibnd=1,nn
          if(ibnd<=max_i) read(iunu) ene_dft_xc_off(1:max_i,ibnd)
       enddo
       close(iunu)
    endif
    call mp_bcast(ene_dft_xc_off, ionode_id, world_comm)


    return
  end SUBROUTINE read_data_pw_dft_xc_off



   SUBROUTINE read_data_pw_upper_states(us,prefix)
!this subroutine reads in the upper REDUCED states

    USE kinds,                ONLY : DP
    USE io_global,            ONLY : stdout, ionode, ionode_id
    USE basic_structures
    USE mp,                   ONLY : mp_bcast
    USE mp_world,             ONLY : world_comm
   USE io_files,  ONLY :  tmp_dir

    implicit none
    INTEGER, EXTERNAL :: find_free_unit
    TYPE(upper_states) :: us!structure to be read and initialized  
    CHARACTER(LEN=256) ::  prefix!to designate the PW files        

    INTEGER :: iunu
    INTEGER :: ii

    if(ionode) then
       iunu = find_free_unit()
       open( unit=iunu, file=trim(tmp_dir)//trim(prefix)//'.upper', status='old',form='unformatted')
       read(iunu) us%nums_tot
       read(iunu) us%nums
       read(iunu) us%nums_occ
       read(iunu) us%nums_reduced
    endif

    call mp_bcast(us%nums_tot, ionode_id, world_comm)
    call mp_bcast(us%nums, ionode_id, world_comm)
    call mp_bcast(us%nums_occ, ionode_id, world_comm)
    call mp_bcast(us%nums_reduced, ionode_id, world_comm)
    
    allocate(us%ene(us%nums_reduced))

    if(ionode) then
       do ii=1,us%nums_reduced
          read(iunu) us%ene(ii)
       enddo
       close(iunu)
    endif
    call mp_bcast(us%ene(:),ionode_id, world_comm)

    return
  END SUBROUTINE read_data_pw_upper_states


SUBROUTINE read_data_pw_vt_mat_lanczos(vtl, ii, prefix, l_pola, ispin)
!this subroutine reads the terms V^v_{v,l}=<Pc w_v(r)w^P_i(r)|z^v_l> from disk
    USE kinds,                ONLY : DP
    USE io_global,            ONLY : stdout, ionode, ionode_id
    USE basic_structures,     ONLY : vt_mat_lanczos,free_memory,initialize_memory
    USE mp,                   ONLY : mp_bcast, mp_barrier, mp_sum
    USE mp_world,             ONLY : world_comm
   USE io_files,  ONLY :  tmp_dir

    implicit none
    
    INTEGER, EXTERNAL :: find_free_unit
    TYPE(vt_mat_lanczos) :: vtl!the structure to be read
    CHARACTER(LEN=256) ::  prefix!to designate the PW files
    LOGICAL ::  l_pola !if true reads the terms for the polarization, otherwise for the self-energy
    INTEGER :: ii!state to be read
    INTEGER, INTENT(in) :: ispin!spin channel


    CHARACTER(4) :: nfile
    INTEGER :: iuntmat, il
    INTEGER, PARAMETER ::  offset=0!ATTENZIONE RESTART it should be 0 normalwise

         
    call initialize_memory(vtl)
    call free_memory(vtl)

    vtl%ii=ii
    write(nfile,'(4i1)') &
       & vtl%ii/1000,mod(vtl%ii,1000)/100,mod(vtl%ii,100)/10,mod(vtl%ii,10)

    if(ionode) then
       iuntmat=find_free_unit()
       if(ispin==1) then
          if(l_pola) then
             open( unit= iuntmat, file=trim(tmp_dir)//trim(prefix)//'.p_mat_lanczos'//nfile, status='old',form='unformatted')
          else
             open( unit= iuntmat, file=trim(tmp_dir)//trim(prefix)//'.s_mat_lanczos'//nfile, status='old',form='unformatted')
          endif
       else
          if(l_pola) then
             open( unit= iuntmat, file=trim(tmp_dir)//trim(prefix)//'.p_mat_lanczos2'//nfile, status='old',form='unformatted')
          else
             open( unit= iuntmat, file=trim(tmp_dir)//trim(prefix)//'.s_mat_lanczos2'//nfile, status='old',form='unformatted')
          endif
       endif
       
       read(iuntmat) vtl%ii
       read(iuntmat) vtl%nums_occ
       read(iuntmat) vtl%numpw
       read(iuntmat) vtl%numl
       vtl%numl=vtl%numl-offset
    endif


    call mp_bcast(vtl%nums_occ,ionode_id, world_comm)
    call mp_bcast(vtl%numpw,ionode_id, world_comm)
    call mp_bcast(vtl%numl,ionode_id, world_comm)


    allocate(vtl%vt_mat(vtl%numpw,vtl%numl))
    if(ionode) then
       do il=1,offset
          read(iuntmat) vtl%vt_mat(1:vtl%numpw,1)
       enddo
    endif
    do il=offset+1,vtl%numl+offset
       !call mp_barrier( world_comm )
       if(ionode) then
          read(iuntmat) vtl%vt_mat(1:vtl%numpw,il-offset)
       else
          vtl%vt_mat(1:vtl%numpw,il-offset)=0.d0
       endif
       !call mp_bcast(vtl%vt_mat(:,il),ionode_id, world_comm)
       !call mp_sum(vtl%vt_mat(1:vtl%numpw,il))
    enddo
    call mp_bcast(vtl%vt_mat,ionode_id, world_comm)
    if(ionode) close(iuntmat)


    return
  END SUBROUTINE read_data_pw_vt_mat_lanczos


  SUBROUTINE read_data_pw_mat_lanczos_full(fl, ii, prefix)
!this subroutine read the full relativistic overlaps

    USE kinds,                ONLY : DP
    USE io_global,            ONLY : stdout, ionode, ionode_id
    USE basic_structures,     ONLY : mat_lanczos_full,free_memory,initialize_memory
    USE mp,                   ONLY : mp_bcast, mp_barrier, mp_sum
    USE mp_world,             ONLY : world_comm
   USE io_files,  ONLY :  tmp_dir

    implicit none

    INTEGER, EXTERNAL :: find_free_unit
    TYPE(mat_lanczos_full) :: fl!the structure to be read
    CHARACTER(LEN=256) ::  prefix!to designate the PW files
    INTEGER :: ii!state to be read


    CHARACTER(4) :: nfile
    INTEGER :: iun, iw,idumm
  
    fl%ii=ii
    write(nfile,'(4i1)') &
       &fl%ii/1000,mod(fl%ii,1000)/100,mod(fl%ii,100)/10,mod(fl%ii,10)

    call initialize_memory(fl)

    if(ionode) then
       iun=find_free_unit()
       open( unit= iun, file=trim(tmp_dir)//trim(prefix)//'.s_mat_lanczos_full'//nfile, status='old',form='unformatted')
       read(iun) idumm
       read(iun) fl%numpw
       read(iun) fl%nums
    endif
    call mp_bcast(fl%numpw, ionode_id, world_comm)
    call mp_bcast(fl%nums, ionode_id, world_comm)
    allocate(fl%f_mat(fl%numpw,fl%nums,2))
    if(ionode) then
       do iw=1,fl%nums
          read(iun) fl%f_mat(1:fl%numpw,iw,1)
       enddo
       do iw=1,fl%nums
          read(iun) fl%f_mat(1:fl%numpw,iw,2)
       enddo
       close(iun)
    endif
    call mp_bcast(fl%f_mat, ionode_id, world_comm)

    return
  END SUBROUTINE read_data_pw_mat_lanczos_full





SUBROUTINE read_data_pw_tt_mat_lanczos(ttl, ii, prefix, l_pola,ispin)
!this subroutine reads the termsT^v_{i,j}=<z^v_i|t_j> from disk
    USE kinds,                ONLY : DP
    USE io_global,            ONLY : stdout, ionode, ionode_id
    USE basic_structures,     ONLY : tt_mat_lanczos,free_memory,initialize_memory
    USE mp,                   ONLY : mp_bcast, mp_barrier, mp_sum
    USE mp_world,             ONLY : world_comm
   USE io_files,  ONLY :  tmp_dir

    implicit none
    INTEGER, EXTERNAL :: find_free_unit

    TYPE(tt_mat_lanczos) :: ttl!the structure to be read
    CHARACTER(LEN=256) ::  prefix!to designate the PW files
    LOGICAL ::  l_pola !if true reads the terms for the polarization, otherwise for the self-energy
    INTEGER :: ii!state to be read
    INTEGER, INTENT(in) :: ispin!spin channel

    CHARACTER(4) :: nfile
    INTEGER :: iuntmat, il


    call initialize_memory(ttl)
    call free_memory(ttl)

    ttl%ii=ii
    write(nfile,'(4i1)') &
       & ttl%ii/1000,mod(ttl%ii,1000)/100,mod(ttl%ii,100)/10,mod(ttl%ii,10)

    if(ionode) then
       iuntmat=find_free_unit()
       if(ispin==1) then
          if(l_pola) then
             open( unit= iuntmat, file=trim(tmp_dir)//trim(prefix)//'.pt_mat_lanczos'//nfile, status='old',form='unformatted')
          else
             open( unit= iuntmat, file=trim(tmp_dir)//trim(prefix)//'.st_mat_lanczos'//nfile, status='old',form='unformatted')
          endif
       else
          if(l_pola) then
             open( unit= iuntmat, file=trim(tmp_dir)//trim(prefix)//'.pt_mat_lanczos2'//nfile, status='old',form='unformatted')
          else
             open( unit= iuntmat, file=trim(tmp_dir)//trim(prefix)//'.st_mat_lanczos2'//nfile, status='old',form='unformatted')
          endif
       endif
       read(iuntmat) ttl%numt
       read(iuntmat) ttl%numl
       read(iuntmat) ttl%ii
    endif

    call mp_bcast(ttl%numt,ionode_id, world_comm)
    call mp_bcast(ttl%numl,ionode_id, world_comm)
  

    allocate(ttl%tt_mat(ttl%numt,ttl%numl))
    do il=1,ttl%numl
       !call mp_barrier
       if(ionode) then
          read(iuntmat) ttl%tt_mat(1:ttl%numt,il)
       else
           ttl%tt_mat(1:ttl%numt,il)=0.d0
        endif
       !call mp_bcast(ttl%tt_mat(:,il),ionode_id, world_comm)
        !call mp_sum( ttl%tt_mat(1:ttl%numt,il))
    enddo
    call mp_bcast(ttl%tt_mat,ionode_id, world_comm)
    if(ionode) close(iuntmat)


    return
  END SUBROUTINE read_data_pw_tt_mat_lanczos

 
SUBROUTINE read_data_pw_lanczos_chain(lc, ii, prefix, l_pola,ispin)
!this subroutine reads the lanczos chain descriptor from disk                                                               
!the date are distributed over the processors
    USE kinds,                ONLY : DP
    USE io_global,            ONLY : stdout, ionode, ionode_id
    USE basic_structures,     ONLY : lanczos_chain,free_memory,initialize_memory
    USE mp,                   ONLY : mp_bcast, mp_barrier, mp_sum
    USE mp_world,             ONLY : nproc,mpime, world_comm
   USE io_files,  ONLY :  tmp_dir

    implicit none
    INTEGER, EXTERNAL :: find_free_unit

    TYPE(lanczos_chain) :: lc!the structure to be read  
    CHARACTER(LEN=256) ::  prefix!to designate the PW files 
    LOGICAL ::  l_pola !if true reads the terms for the polarization, otherwise for the self-energy 
    INTEGER :: ii!state to be read , only for self-energy
    INTEGER, INTENT(in) :: ispin!spin multiplicity

    CHARACTER(4) :: nfile
    INTEGER :: iunlc, is,it
    INTEGER :: l_blk,nbegin,nend
    REAL(kind=DP), ALLOCATABLE :: tmp_mat(:)
    

    call initialize_memory(lc)
    call free_memory(lc)

    lc%ii=ii
    write(nfile,'(4i1)') &
       & lc%ii/1000,mod(lc%ii,1000)/100,mod(lc%ii,100)/10,mod(lc%ii,10)

    if(ionode) then
       iunlc=find_free_unit()
       if(ispin==1) then
          if(l_pola) then
             open( unit= iunlc, file=trim(tmp_dir)//trim(prefix)//'.p_iter_lanczos', status='old',form='unformatted')
          else
             open( unit= iunlc, file=trim(tmp_dir)//trim(prefix)//'.s_iter_lanczos'//'_'//nfile, status='old',form='unformatted')
          endif
       else
          if(l_pola) then
             open( unit= iunlc, file=trim(tmp_dir)//trim(prefix)//'.p_iter_lanczos2', status='old',form='unformatted')
          else
             open( unit= iunlc, file=trim(tmp_dir)//trim(prefix)//'.s_iter_lanczos2'//'_'//nfile, status='old',form='unformatted')
          endif
       endif
       read(iunlc) lc%numt
       read(iunlc) lc%ii
       read(iunlc) lc%num_steps
    endif

    write(*,*) lc%numt,  lc%ii,lc%num_steps
    call mp_bcast(lc%numt,ionode_id, world_comm)
    call mp_bcast(lc%num_steps,ionode_id, world_comm)
    
    l_blk= (lc%numt)/nproc   
    if(l_blk*nproc < (lc%numt)) l_blk = l_blk+1
    nbegin=mpime*l_blk+1
    nend=nbegin+l_blk-1
    allocate(tmp_mat(lc%numt))

    
    allocate(lc%o_mat(lc%numt,lc%num_steps,l_blk))
    allocate(lc%d(lc%num_steps,lc%numt))
    allocate(lc%f(lc%num_steps,lc%numt))

    do is=1,lc%num_steps
       do it=1,lc%numt
          tmp_mat(1:lc%numt)=0.d0
          if(ionode) read(iunlc)  tmp_mat(1:lc%numt)
          call mp_sum(tmp_mat(:),world_comm)!this should be faster than mp_bcat
          if(it >= nbegin .and. it <= nend) then
             lc%o_mat(1:lc%numt,is,it-nbegin+1)= tmp_mat(1:lc%numt)
          endif
         ! if(ionode) read(iunlc) lc%o_mat(1:lc%numt,is,it)
         ! call mp_barrier
         ! call mp_bcast(lc%o_mat(1:lc%numt,is,it),ionode_id, world_comm)
       enddo
    enddo

    do it=1,lc%numt
       if(ionode) read(iunlc) lc%d(1:lc%num_steps,it)
       call mp_barrier( world_comm )
       call mp_bcast(lc%d(1:lc%num_steps,it),ionode_id, world_comm)
    enddo

      do it=1,lc%numt
       if(ionode) read(iunlc) lc%f(1:lc%num_steps,it)
       call mp_barrier( world_comm )
       call mp_bcast(lc%f(1:lc%num_steps,it),ionode_id, world_comm)
    enddo

    if(ionode) close(iunlc)

    deallocate(tmp_mat)

    return

  END SUBROUTINE read_data_pw_lanczos_chain


SUBROUTINE read_data_pw_vt_mat_lanczos_single(vtl, ii, prefix, l_pola)
!this subroutine reads the terms V^v_{v,l}=<Pc w_v(r)w^P_i(r)|z^v_l> from disk
!single processor version
    USE kinds,                ONLY : DP
    USE basic_structures,     ONLY : vt_mat_lanczos,free_memory,initialize_memory
   USE io_files,  ONLY :  tmp_dir


    implicit none
    INTEGER, EXTERNAL :: find_free_unit
    
    TYPE(vt_mat_lanczos) :: vtl!the structure to be read
    CHARACTER(LEN=256) ::  prefix!to designate the PW files
    LOGICAL ::  l_pola !if true reads the terms for the polarization, otherwise for the self-energy
    INTEGER :: ii!state to be read


    CHARACTER(4) :: nfile
    INTEGER :: iuntmat, il
    INTEGER, PARAMETER ::  offset=0!ATTENZIONE RESTART it should be 0 normalwise
    
         
    call initialize_memory(vtl)
    call free_memory(vtl)

    vtl%ii=ii
    write(nfile,'(4i1)') &
       & vtl%ii/1000,mod(vtl%ii,1000)/100,mod(vtl%ii,100)/10,mod(vtl%ii,10)


    iuntmat=find_free_unit()
    if(l_pola) then
       open( unit= iuntmat, file=trim(tmp_dir)//trim(prefix)//'.p_mat_lanczos'//nfile, status='old',form='unformatted')
    else
       open( unit= iuntmat, file=trim(tmp_dir)//trim(prefix)//'.s_mat_lanczos'//nfile, status='old',form='unformatted')
    endif
       
    read(iuntmat) vtl%ii
    read(iuntmat) vtl%nums_occ
    read(iuntmat) vtl%numpw
    read(iuntmat) vtl%numl
    vtl%numl=vtl%numl-offset

    
    allocate(vtl%vt_mat(vtl%numpw,vtl%numl))
    do il=1,offset
       read(iuntmat) vtl%vt_mat(1:vtl%numpw,1)
    enddo
    do il=1+offset,vtl%numl+offset
       read(iuntmat) vtl%vt_mat(1:vtl%numpw,il-offset)
    enddo
    
             
    close(iuntmat)


    return
  END SUBROUTINE read_data_pw_vt_mat_lanczos_single


SUBROUTINE read_data_pw_tt_mat_lanczos_single(ttl, ii, prefix, l_pola)
!this subroutine reads the termsT^v_{i,j}=<z^v_i|t_j> from disk
!single processor version
    USE kinds,                ONLY : DP
    USE basic_structures,     ONLY : tt_mat_lanczos,free_memory,initialize_memory
    USE mp,                   ONLY : mp_bcast, mp_barrier
    USE mp_world,             ONLY : world_comm
   USE io_files,  ONLY :  tmp_dir

    implicit none

    INTEGER, EXTERNAL :: find_free_unit
    TYPE(tt_mat_lanczos) :: ttl!the structure to be read
    CHARACTER(LEN=256) ::  prefix!to designate the PW files
    LOGICAL ::  l_pola !if true reads the terms for the polarization, otherwise for the self-energy
    INTEGER :: ii!state to be read


    CHARACTER(4) :: nfile
    INTEGER :: iuntmat, il


    call initialize_memory(ttl)
    call free_memory(ttl)

    ttl%ii=ii
    write(nfile,'(4i1)') &
       & ttl%ii/1000,mod(ttl%ii,1000)/100,mod(ttl%ii,100)/10,mod(ttl%ii,10)


    iuntmat=find_free_unit()
    if(l_pola) then
       open( unit= iuntmat, file=trim(tmp_dir)//trim(prefix)//'.pt_mat_lanczos'//nfile, status='old',form='unformatted')
    else
       open( unit= iuntmat, file=trim(tmp_dir)//trim(prefix)//'.st_mat_lanczos'//nfile, status='old',form='unformatted')
    endif
    
    read(iuntmat) ttl%numt
    read(iuntmat) ttl%numl
    read(iuntmat) ttl%ii


    
  

    allocate(ttl%tt_mat(ttl%numt,ttl%numl))
    do il=1,ttl%numl
       read(iuntmat) ttl%tt_mat(1:ttl%numt,il)
    enddo

    close(iuntmat)


    return
  END SUBROUTINE read_data_pw_tt_mat_lanczos_single


  SUBROUTINE read_data_pw_full_prods(fp,prefix)
!this subroutine read the full relativistic overlaps                                                                                            

    USE kinds,                ONLY : DP
    USE io_global,            ONLY : stdout, ionode, ionode_id
    USE basic_structures,     ONLY : full_prods,free_memory,initialize_memory
    USE mp,                   ONLY : mp_bcast, mp_barrier, mp_sum
    USE mp_world,             ONLY : world_comm
   USE io_files,  ONLY :  tmp_dir

    implicit none

    INTEGER, EXTERNAL :: find_free_unit
    TYPE(full_prods) :: fp!the structure to be read 
    CHARACTER(LEN=256) ::  prefix!to designate the PW files

    INTEGER :: iun, is,ii,ipol

    iun=find_free_unit()
    if(ionode) then
       open( unit= iun, file=trim(tmp_dir)//trim(prefix)//'.prod_full', status='old',form='unformatted')
       read(iun) fp%nums
       read(iun) fp%nbnd
       read(iun) fp%numpw
       read(iun) fp%numv
    endif
    call mp_bcast(fp%nums,ionode_id, world_comm)
    call mp_bcast(fp%nbnd, ionode_id, world_comm)
    call mp_bcast(fp%numpw, ionode_id, world_comm)
    call mp_bcast(fp%numv, ionode_id, world_comm)
    allocate(fp%ene_ks(fp%nbnd))
    allocate(fp%gmat(fp%numpw,2,fp%nbnd,fp%nums))
    if(ionode) then
       read(iun) fp%ene_ks(1:fp%nbnd)
       do is=1,fp%nums
          do ii=1,fp%nbnd
             do ipol=1,2
                read(iun) fp%gmat(1:fp%numpw,ipol,ii,is)
             enddo
          enddo
       enddo
       close(iun)
    endif

    call mp_bcast(fp%ene_ks, ionode_id, world_comm)
    call mp_bcast(fp%gmat, ionode_id, world_comm)
    
    return
  END SUBROUTINE read_data_pw_full_prods


SUBROUTINE read_data_pw_partial_occ(po, prefix, ispin)

  USE kinds,                ONLY : DP
  USE basic_structures,     ONLY : partial_occ,free_memory,initialize_memory
  USE mp,                   ONLY : mp_bcast, mp_barrier
  USE mp_world,             ONLY : world_comm
  USE io_global,            ONLY : ionode, ionode_id
   USE io_files,  ONLY :  tmp_dir

  implicit none
  INTEGER, EXTERNAL :: find_free_unit
  TYPE(partial_occ), INTENT(out) :: po!the structure to be read                                                                 
  CHARACTER(LEN=256),INTENT(in) ::  prefix!to designate the PW files                                                
  INTEGER, INTENT(in) :: ispin!spin channel
  
  INTEGER :: iun, iv,jv

  call free_memory(po)
  if(ionode) then
     iun=find_free_unit()
     if(ispin==1) then
        open( unit= iun, file=trim(tmp_dir)//trim(prefix)//'.occ_mat', status='old',form='unformatted')
     else
        open( unit= iun, file=trim(tmp_dir)//trim(prefix)//'.occ_mat2', status='old',form='unformatted')
     endif
     read(iun) po%nums_occ_min
     read(iun) po%nums_occ
     read(iun) po%numpw
  endif
  call mp_bcast(po%nums_occ_min,ionode_id, world_comm)
  call mp_bcast(po%nums_occ, ionode_id, world_comm)
  call mp_bcast(po%numpw, ionode_id, world_comm)
  allocate(po%f_occ(po%nums_occ))
  if(ionode) read(iun) po%f_occ(1:po%nums_occ)
  call mp_bcast(po%f_occ, ionode_id, world_comm)
  allocate(po%ppp_mat(po%numpw,po%nums_occ,po%nums_occ_min+1:po%nums_occ))
  do iv=po%nums_occ_min+1,po%nums_occ
     do jv=1,po%nums_occ
        if(ionode) read(iun) po%ppp_mat(1:po%numpw,jv,iv)
        call mp_bcast( po%ppp_mat(1:po%numpw,jv,iv),ionode_id, world_comm)
     enddo
  enddo
  if(ionode) close(iun)

  return
  
END SUBROUTINE read_data_pw_partial_occ


SUBROUTINE read_data_pw_semicore(sc, prefix, ispin)
!NOT_TO_BE_INCLUDED_START
  USE kinds,                ONLY : DP
  USE basic_structures,     ONLY : semicore,free_memory,initialize_memory
  USE mp,                   ONLY : mp_bcast, mp_barrier
  USE mp_world,             ONLY : world_comm
  USE io_global,            ONLY : ionode, ionode_id
   USE io_files,  ONLY :  tmp_dir

  implicit none
  INTEGER, EXTERNAL :: find_free_unit
  TYPE(semicore), INTENT(out) :: sc!the structure to be read
  CHARACTER(LEN=256),INTENT(in) ::  prefix!to designate the PW files
  INTEGER, INTENT(in) :: ispin!spin channel

  INTEGER :: iun, iw,ii
  REAL(kind=DP), ALLOCATABLE :: tmp_prod(:)

  if(ionode) then
     iun=find_free_unit()
     if(ispin==1) then
        open( unit= iun, file=trim(tmp_dir)//trim(prefix)//'.sc_gvphi', status='old',form='unformatted')
     else
        open( unit= iun, file=trim(tmp_dir)//trim(prefix)//'.sc_gvphi2', status='old',form='unformatted')
     endif
     read(iun) sc%n_semicore
  endif
  call mp_bcast(sc%n_semicore, ionode_id, world_comm)
  allocate(sc%en_sc(sc%n_semicore))
  if(ionode) then
     read(iun) sc%en_sc(1:sc%n_semicore)
     read(iun) sc%nums
     read(iun) sc%numpw
  endif
  call mp_bcast(sc%en_sc,ionode_id, world_comm)
  call mp_bcast(sc%nums, ionode_id, world_comm)
  call mp_bcast(sc%numpw, ionode_id, world_comm)

  allocate(sc%ppw_mat(sc%numpw,sc%n_semicore,sc%nums))

  allocate(tmp_prod(sc%n_semicore))
  if(ionode) then
     do iw=1,sc%numpw
        do ii=1,sc%nums
           read(iun) tmp_prod(1:sc%n_semicore)
           sc%ppw_mat(iw,1:sc%n_semicore,ii)= tmp_prod(1:sc%n_semicore)
        enddo
     enddo
  endif
  call mp_bcast(sc%ppw_mat, ionode_id, world_comm)

  deallocate(tmp_prod)
  if(ionode) close(iun)

  return

!NOT_TO_BE_INCLUDED_END
END SUBROUTINE read_data_pw_semicore


SUBROUTINE read_data_pw_contour(ct,prefix,ispin,istate)
!NOT_TO_BE_INCLUDED_START
!this subroutines reads the <psi_i|s_\alpha> overlaps


  USE kinds,                ONLY : DP
  USE basic_structures,     ONLY : contour_terms,free_memory,initialize_memory
  USE mp,                   ONLY : mp_bcast, mp_barrier
  USE mp_world,             ONLY : world_comm
  USE io_global,            ONLY : ionode, ionode_id
  USE io_files,  ONLY :  tmp_dir

  implicit none
  INTEGER, EXTERNAL :: find_free_unit
  TYPE(contour_terms), INTENT(out) :: ct!the structure to be read
  CHARACTER(LEN=256),INTENT(in) ::  prefix!to designate the PW files    
  INTEGER, INTENT(in) :: ispin!spin channel
  INTEGER, INTENT(in) :: istate!!KS states relative to global s vectors for big_system option

  INTEGER :: iun, iw,ii
  CHARACTER(4) :: nfile

  if(ionode) then
     iun=find_free_unit()
      write(nfile,'(4i1)') istate/1000,mod(istate,1000)/100,mod(istate,100)/10,mod(istate,10)
     if(ispin==1) then
        open( unit= iun, file=trim(tmp_dir)//trim(prefix)//'.s_contour'//nfile, status='old',form='unformatted')
     else
        open( unit= iun, file=trim(tmp_dir)//trim(prefix)//'.s_contour'//nfile, status='old',form='unformatted')
     endif
     read(iun) ct%nums
     read(iun) ct%numt
  endif
  call mp_bcast(ct%nums, ionode_id, world_comm)
  call mp_bcast(ct%numt, ionode_id, world_comm)
  allocate(ct%cmat(ct%numt,ct%nums))
  if(ionode) then
     do ii=1,ct%nums
        read(iun) ct%cmat(1:ct%numt,ii)
     enddo
     close(iun)
  endif
   call mp_bcast(ct%cmat, ionode_id, world_comm)

  return
!NOT_TO_BE_INCLUDED_END
END SUBROUTINE read_data_pw_contour
