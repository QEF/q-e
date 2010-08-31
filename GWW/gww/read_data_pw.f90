!P.Umari  GWW code
! Modified by G. Stenuit
!
!these subroutines read in the data from PW calculations



SUBROUTINE read_data_pw_u(wu,prefix)
!this subroutine reads in the energies and the inversewannier transformation matrix

    USE kinds,                ONLY : DP
    USE io_global,            ONLY : stdout, ionode, ionode_id
    USE io_files,             ONLY : find_free_unit
    USE basic_structures
    USE mp,                   ONLY : mp_bcast, mp_barrier

    implicit none

    TYPE(wannier_u) :: wu!structure to be read and initialized
    CHARACTER(LEN=256) ::  prefix!to designate the PW files

    LOGICAL :: lex
    INTEGER :: iunu, nn
    INTEGER :: iw

    if(ionode) then
       iunu = find_free_unit()
       open( unit=iunu, file=trim(prefix)//'.wannier', status='old',form='unformatted')

!read in basis length

       read(iunu) wu%nums
       read(iunu) wu%nums_occ
    endif

    call mp_bcast(wu%nums, ionode_id)
    call mp_bcast(wu%nums_occ, ionode_id)

!allocate arrays
    allocate(wu%ene(wu%nums))
    allocate(wu%ene_xc(wu%nums))
    allocate(wu%ene_u(wu%nums))
    allocate(wu%ene_lda_h(wu%nums))
    allocate(wu%umat(wu%nums,wu%nums))

    if(ionode) then
!read in energies
       read(iunu) wu%ene(1:wu%nums)
!read in DFT exchange and correlation energies

       read(iunu) wu%ene_xc(1:wu%nums)
       read(iunu) wu%ene_lda_h(1:wu%nums)
!read in transformation matrix
       do iw=1,wu%nums
          read(iunu) wu%umat(1:wu%nums,iw)
       enddo
    endif

    call mp_bcast(wu%ene(:), ionode_id)
    call mp_bcast(wu%ene_xc(:), ionode_id)
    call mp_bcast(wu%ene_lda_h(:), ionode_id)
    do iw=1,wu%nums
       call mp_barrier
       call mp_bcast(wu%umat(:,iw), ionode_id)
    enddo

    if(ionode) close(iunu)

!!!!!!!!! add the contribution due to the U
    if(ionode) then
       !
       INQUIRE ( file=trim(prefix)//'.hubbard_u', EXIST=lex )
       write(stdout,*) 'hubbard_u files exist ? lex=', lex
       call flush_unit(stdout)
       !
       if ( lex ) then
         !
         iunu = find_free_unit()
         open( unit=iunu, file=trim(prefix)//'.hubbard_u', status='old',form='unformatted')
         read(iunu) nn
         do iw=1,wu%nums
            read(iunu) wu%ene_u(iw)
         enddo
         !
       else
         !
         do iw=1,wu%nums
            wu%ene_u(iw)=0.0
         enddo
         !
       endif
    endif
    !
    call mp_bcast(wu%ene_u(:), ionode_id)
    !
    if(ionode) close(iunu)
    !
    ! MODIFY the xc and add the U contribution:
    do iw=1,wu%nums
       write(stdout,*) 'OLD wu%ene_xc(',iw,')=', wu%ene_xc(iw)
       write(stdout,*) 'wu%ene_u(',iw,')=', wu%ene_u(iw)
    enddo
    call flush_unit(stdout)
    !
    wu%ene_xc=wu%ene_xc+wu%ene_u
    !
    do iw=1,wu%nums
       write(stdout,*) 'NEW wu%ene_xc(',iw,')=', wu%ene_xc(iw)
    enddo
    call flush_unit(stdout)
    !
    return
    !
END SUBROUTINE read_data_pw_u


SUBROUTINE  read_data_pw_v(vp,prefix,debug,ort,l_zero)
!read from file and initialize coulomb potential on the basis of products of wanniers

    USE kinds,                ONLY : DP
    USE io_global,            ONLY : stdout, ionode, ionode_id
    USE io_files,             ONLY : find_free_unit
    USE constants,            ONLY : eps8
    USE basic_structures
    USE mp,                   ONLY : mp_bcast,mp_barrier

    implicit none

    TYPE(v_pot) :: vp!structure to be read and initialized
    CHARACTER(LEN=256) ::  prefix!to designate the PW files
    LOGICAL :: debug!if true check for simmetry
    INTEGER :: ort!if ort==0, open non orthogonal file, if ort == 1 open orthogonal file
                  !if ort==2 open non orthogonal symmetric file
    LOGICAL :: l_zero!if true open file with head put to zero

    INTEGER :: iunv
    INTEGER :: iw,jw

    write(stdout,*) 'SUBROUTINE read_data_pw_v'
    call flush_unit(stdout)

    if(ionode) then
       iunv=find_free_unit()
       if(ort==1) then
          open( unit=iunv, file=trim(prefix)//'.vpot', status='old',form='unformatted')
       else if (ort==0) then
          if(.not.l_zero) then
             open( unit=iunv, file=trim(prefix)//'.vpot_no', status='old',form='unformatted')
          else
             open( unit=iunv, file=trim(prefix)//'.vpot_no_zero', status='old',form='unformatted')
          endif
       else if (ort==2) then
          if(.not.l_zero) then
             open( unit=iunv, file=trim(prefix)//'.vpot_no_sym', status='old',form='unformatted')
          else
             open( unit=iunv, file=trim(prefix)//'.vpot_no_sym_zero', status='old',form='unformatted')
          endif
       endif
!read in basis length
       read(iunv) vp%numpw

    endif

    write(stdout,*) 'After reading the .vpot_ files: vp%numpw=', vp%numpw
    call flush_unit(stdout)

    call mp_bcast(vp%numpw, ionode_id)

    write(stdout,*) 'After the first mp_bcast'
    call flush_unit(stdout)

!allocate array
    allocate(vp%vmat(vp%numpw,vp%numpw))

!read in potential matrix
    if(ionode) then
       do iw=1,vp%numpw
          read(iunv) vp%vmat(1:vp%numpw,iw)
       enddo
    endif

    write(stdout,*) 'After reading the potential matrix'
    call flush_unit(stdout)

    do iw=1,vp%numpw
       call mp_barrier
       call mp_bcast(vp%vmat(:,iw), ionode_id)
    enddo

    write(stdout,*) 'After the seconf mp_bcast'
    call flush_unit(stdout)

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

    write(stdout,*) 'END SUBROUTINE read_data_pw_v'
    call flush_unit(stdout)

    return
END SUBROUTINE read_data_pw_v


SUBROUTINE read_data_pw_q(qm,prefix,l_v_products)
!this subroutine reads in and allocate the arrays for the
!description of overlaps of (orthonormalized) products of wanniers
!with products of wannier


    USE kinds,                ONLY : DP
    USE io_global,            ONLY : stdout, ionode, ionode_id
    USE io_files,             ONLY : find_free_unit
    USE basic_structures
    USE mp,                   ONLY : mp_bcast, mp_barrier

    implicit none

    TYPE(q_mat) :: qm!structure to be read and initialized
    CHARACTER(LEN=256) ::  prefix!to designate the PW files
    LOGICAL :: l_v_products!if true read the wp_v file for the products \tilde{w}_i(r)\tilde{w}_j(r)v(r,r')w^P_red_k(r')

    INTEGER :: iunq
    INTEGER :: iw


    if(ionode) then
       iunq=find_free_unit()
       if(.not.l_v_products) then
          open( unit=iunq, file=trim(prefix)//'.wp', status='old',form='unformatted')
       else
          open( unit=iunq, file=trim(prefix)//'.wp_v', status='old',form='unformatted')
       endif

!read in basis length
       read(iunq) qm%numpw
    endif

    call mp_bcast(qm%numpw, ionode_id)

! allocate array of descriptors
    allocate (qm%wp(qm%numpw))


    do iw=1,qm%numpw

       if(ionode)    read(iunq) qm%wp(iw)%numij
       call mp_bcast(qm%wp(iw)%numij, ionode_id)

!for each descriptor allocates arrays
       allocate(qm%wp(iw)%ij(2,qm%wp(iw)%numij))
       allocate(qm%wp(iw)%o(qm%wp(iw)%numij))

!read data
       if(ionode) then
          read(iunq) qm%wp(iw)%ij(1,1:qm%wp(iw)%numij)
          read(iunq) qm%wp(iw)%ij(2,1:qm%wp(iw)%numij)
          read(iunq) qm%wp(iw)%o(1:qm%wp(iw)%numij)
       end if
       call mp_bcast(qm%wp(iw)%ij(:,:), ionode_id)
       call mp_bcast(qm%wp(iw)%o(:), ionode_id)
    enddo


    qm%is_parallel=.false.
    qm%numpw_para=qm%numpw
    qm%first_para=1

    if(ionode) close(iunq)

END SUBROUTINE read_data_pw_q


SUBROUTINE read_data_pw_ortho_polaw(op,prefix)
!this subroutine reads in and allocate the arrays for the
!description of orthonormalization matrix of wanniers products


    USE kinds,                ONLY : DP
    USE io_global,            ONLY : stdout, ionode, ionode_id
    USE io_files,             ONLY : find_free_unit
    USE basic_structures,     ONLY : ortho_polaw, free_memory
    USE mp,                   ONLY : mp_bcast, mp_barrier

    implicit none

    TYPE(ortho_polaw) :: op!structure to be read and initialized
    CHARACTER(LEN=256) ::  prefix!to designate the PW files

    INTEGER :: iunq
    INTEGER :: iw

!    call free_memory(op)

    if(ionode) then
       iunq=find_free_unit()
       open( unit=iunq, file=trim(prefix)//'.orthonorm', status='old',form='unformatted')

!read in basis length
       read(iunq) op%numpw
    endif
    call mp_bcast(op%numpw, ionode_id)
    allocate(op%on_mat(op%numpw,op%numpw))

    if(ionode) then
       do iw=1,op%numpw
          read(iunq) op%on_mat(1:op%numpw,iw)
       enddo
    end if
    do iw=1,op%numpw
       call mp_barrier
       call mp_bcast( op%on_mat(:,iw), ionode_id)
    enddo

    op%inverse=.false.

    if(ionode) close(iunq)

END SUBROUTINE read_data_pw_ortho_polaw


SUBROUTINE read_data_pw_wp_psi(wp,prefix)
!this subroutine reads in and allocate the arrays for the
!description of products of valence^2 times two wannier products

    USE kinds,                ONLY : DP
    USE io_global,            ONLY : stdout, ionode, ionode_id
    USE io_files,             ONLY : find_free_unit
    USE basic_structures,     ONLY : wp_psi, free_memory
    USE mp,                   ONLY : mp_bcast

    implicit none

    TYPE(wp_psi) :: wp!structure to be read and initialized
    CHARACTER(LEN=256) ::  prefix!to designate the PW files

    INTEGER :: iunq
    INTEGER :: iw,hw, jw


!    call free_memory(wp)

    if(ionode) then
       iunq=find_free_unit()
       open( unit=iunq, file=trim(prefix)//'.wpwp_psi', status='old',form='unformatted')

!read in basis length
       read(iunq) wp%numpw
       read(iunq) wp%nums_psi
    endif
    call mp_bcast(wp%numpw, ionode_id)
    call mp_bcast(wp%nums_psi, ionode_id)
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
        call mp_bcast( wp%wwp(:,:,hw), ionode_id)
     enddo

    if(ionode) close(iunq)

END SUBROUTINE read_data_pw_wp_psi


SUBROUTINE read_data_pw_u_prim(wu,prefix)
!this subroutine reads in the inverse wannier transformation matrix

    USE kinds,                ONLY : DP
    USE io_global,            ONLY : stdout, ionode, ionode_id
    USE io_files,             ONLY : find_free_unit
    USE basic_structures
    USE mp,                   ONLY : mp_bcast, mp_barrier

    implicit none

    TYPE(wannier_u_prim) :: wu!structure to be read and initialized
    CHARACTER(LEN=256) ::  prefix!to designate the PW files

    INTEGER :: iunu
    INTEGER :: iw

    if(ionode) then
       iunu = find_free_unit()
       open( unit=iunu, file=trim(prefix)//'.wannier_prim', status='old',form='unformatted')

!read in basis length
       read(iunu) wu%nums_prim
       read(iunu) wu%nums_occ
       read(iunu) wu%nums
       write(*,*) 'read_data_pw_u_prim',wu%nums_prim,wu%nums_occ,wu%nums
    endif

    call mp_bcast(wu%nums_prim, ionode_id)
    call mp_bcast(wu%nums_occ, ionode_id)
    call mp_bcast(wu%nums, ionode_id)

!allocate arrays
    allocate(wu%umat(wu%nums_prim,wu%nums_prim))

    if(ionode) then
!read in transformation matrix
       do iw=1,wu%nums_prim
          read(iunu) wu%umat(1:wu%nums_prim,iw)
       enddo
    endif

    do iw=1,wu%nums_prim
       call mp_barrier
       call mp_bcast(wu%umat(:,iw), ionode_id)
    enddo

    if(ionode) close(iunu)
    return

END SUBROUTINE read_data_pw_u_prim


SUBROUTINE read_data_pw_v_pot_prim(vp,prefix, l_zero)
!this subroutine reads in the coulombian potential and the overlap index

    USE kinds,                ONLY : DP
    USE io_global,            ONLY : stdout, ionode, ionode_id
    USE io_files,             ONLY : find_free_unit
    USE basic_structures
    USE mp,                   ONLY : mp_bcast, mp_barrier

    implicit none

    TYPE(v_pot_prim) :: vp!structure to be read and initialized
    CHARACTER(LEN=256) ::  prefix!to designate the PW files
    LOGICAL :: l_zero!if true opens file with head of v pu to zero

    INTEGER :: iunu
    INTEGER :: iw

    if(ionode) then
       iunu = find_free_unit()
       if(.not. l_zero) then
          open( unit=iunu, file=trim(prefix)//'.uterms_prim', status='old',form='unformatted')
       else
          open( unit=iunu, file=trim(prefix)//'.uterms_prim_zero', status='old',form='unformatted')
       endif

!read in basis length
       read(iunu) vp%numpw_prim
       read(iunu) vp%numpw
       write(*,*) 'read_data_pw_v_pot_prim', vp%numpw_prim,vp%numpw
    endif

    call mp_bcast(vp%numpw, ionode_id)
    call mp_bcast(vp%numpw_prim, ionode_id)

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
       call mp_barrier
       call mp_bcast(vp%vmat(:,iw), ionode_id)
    enddo

    if(ionode) then
       iunu = find_free_unit()
       open( unit=iunu, file=trim(prefix)//'.ij_prim', status='old',form='unformatted')
       do iw=1,vp%numpw_prim
          read(iunu) vp%ij(1,iw),vp%ij(2,iw)
       enddo
       close(iunu)
    endif


    call mp_bcast(vp%ij(:,:), ionode_id)

    vp%is_parallel=.false.
    vp%numpw_para=vp%numpw
    vp%first_para=1

END SUBROUTINE read_data_pw_v_pot_prim


SUBROUTINE read_data_pw_wp_psi_cutoff_index(wpi,prefix)
!this subroutine reads in and allocate the arrays for the
!indices describing of products of valence^2 times two wannier products
!when a cutoff is applied

    USE kinds,                ONLY : DP
    USE io_global,            ONLY : stdout, ionode, ionode_id
    USE io_files,             ONLY : find_free_unit
    USE basic_structures,     ONLY : wp_psi_cutoff_index, free_memory
    USE mp,                   ONLY : mp_bcast

    implicit none

    TYPE(wp_psi_cutoff_index) :: wpi!structure to be read and initialized
    CHARACTER(LEN=256) ::  prefix!to designate the PW files

    INTEGER :: iuni
    INTEGER :: i

    if(ionode) then
       iuni=find_free_unit()
       open( unit=iuni, file=trim(prefix)//'.wpwp_psi_index', status='old',form='unformatted')
!read in basis length
       read(iuni) wpi%numpw
       read(iuni) wpi%nums_psi
       read(iuni) wpi%numpwpw
    endif

    call mp_bcast(wpi%numpw, ionode_id)
    call mp_bcast(wpi%nums_psi, ionode_id)
    call mp_bcast(wpi%numpwpw, ionode_id)

    allocate(wpi%index(2,wpi%numpwpw))

    if(ionode) then
       do i=1,wpi%numpwpw
          read(iuni) wpi%index(1,i),wpi%index(2,i)
       enddo
       close(iuni)
    endif

    call mp_bcast(wpi%index, ionode_id)

    return
END SUBROUTINE read_data_pw_wp_psi_cutoff_index


SUBROUTINE read_data_pw_wp_psi_cutoff_data(wpi,wp,prefix)
!this subroutine reads in and allocate the arrays for the
!products of valence^2 times two wannier products when a cutoff is applied

    USE kinds,                ONLY : DP
    USE io_global,            ONLY : stdout, ionode, ionode_id
    USE io_files,             ONLY : find_free_unit
    USE basic_structures,     ONLY : wp_psi_cutoff_index, wp_psi_cutoff_data,free_memory
    USE mp,                   ONLY : mp_bcast

    implicit none

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
       open( unit=iund, file=trim(prefix)//'.wpwp_psi', status='old',form='unformatted')
!read in basis length
       do i=1,wp%nums_psi*wp%numpwpw
          read(iund) pos,state,w
          wp%wwp(pos,state)=w
       enddo
       close(iund)
    endif

    do i=1,wp%nums_psi
       call mp_bcast(wp%wwp(:,i), ionode_id)
    enddo

END SUBROUTINE read_data_pw_wp_psi_cutoff_data


SUBROUTINE read_data_pw_exchange(ene_x,max_i,prefix)
!this subroutine reads in the exchange energies

    USE kinds,                ONLY : DP
    USE io_global,            ONLY : stdout, ionode, ionode_id
    USE io_files,             ONLY : find_free_unit
    USE basic_structures
    USE mp,                   ONLY : mp_bcast

    implicit none

    REAL(kind=DP) :: ene_x(max_i)
    INTEGER :: max_i
    CHARACTER(LEN=256) ::  prefix!to designate the PW files

    INTEGER :: iunu
    INTEGER :: ndata
    REAL(kind=DP), ALLOCATABLE :: buf(:)

    if(ionode) then
       iunu = find_free_unit()
       open( unit=iunu, file=trim(prefix)//'.exchange', status='old',form='unformatted')
       read(iunu) ndata
       allocate(buf(ndata))
       read(iunu) buf(1:ndata)
       close(iunu)
       ene_x(1:max_i)=buf(1:max_i)
       deallocate(buf)
    endif
    call mp_bcast(ene_x, ionode_id)

END SUBROUTINE read_data_pw_exchange


SUBROUTINE read_data_pw_head_epsilon(he, prefix, l_wing_epsilon)
!this subroutine reads the data


    USE kinds,                ONLY : DP
    USE io_global,            ONLY : stdout, ionode, ionode_id
    USE io_files,             ONLY : find_free_unit
    USE basic_structures,     ONLY : head_epsilon
    USE mp,                   ONLY : mp_bcast, mp_barrier

    implicit none

    TYPE(head_epsilon) :: he!the head of epsilon to be read
    CHARACTER(LEN=256) ::  prefix!to designate the PW files
    LOGICAL :: l_wing_epsilon!if true read from file also the wing data

    INTEGER :: iun,i, idumm
    REAL(kind=DP) :: rdumm

    if(ionode) then
       iun = find_free_unit()
       open( unit=iun, file=trim(prefix)//'.head', status='old',form='unformatted')
       read(iun) he%n
       read(iun) he%omega
    endif
    call mp_bcast(he%n, ionode_id)
    call mp_bcast(he%omega, ionode_id)
    allocate(he%freqs(he%n+1))
    allocate(he%head(he%n+1))
    if(ionode) then
       read(iun) he%freqs(1:he%n+1)
       read(iun) he%head(1:he%n+1)
       close(iun)
    endif
    call mp_bcast(he%freqs, ionode_id)
    call mp_bcast(he%head, ionode_id)

    if(ionode) then
       iun = find_free_unit()
       open( unit=iun, file=trim(prefix)//'.gzero', status='old',form='unformatted')
       read(iun) he%numpw
    endif
    call mp_bcast(he%numpw, ionode_id)
    allocate(he%gzero(he%numpw))
    if(ionode) then
       do i=1,he%numpw
          read(iun) he%gzero(i)
       enddo
       close(iun)
    endif
    call mp_bcast(he%gzero,ionode_id)

    if(l_wing_epsilon) then
       allocate(he%wing(he%numpw, he%n+1))
       allocate(he%wing_c(he%numpw, he%n+1))
       if(ionode) then
          iun = find_free_unit()
          open( unit=iun, file=trim(prefix)//'.wing', status='old',form='unformatted')
          read(iun) idumm
          read(iun) rdumm
          if(idumm /= he%n) then
             write(stdout,*) 'WING: PROBLEM WITH N'
          endif
          if(rdumm /= he%omega) then
             write(stdout,*) 'WING: PROBLEM WITH OMEGA'
          endif
          read(iun) idumm
          if(idumm /= he%numpw) then
             write(stdout,*) 'WING: PROBLEM WITH NUMPW', idumm, he%numpw
          endif
          do i=1,he%n+1
             read(iun) he%wing(1:he%numpw,i)
          enddo
          do i=1,he%n+1
             read(iun) he%wing_c(1:he%numpw,i)
          enddo

          close(iun)
       endif
       do i=1,he%n+1
          call mp_barrier
          call mp_bcast(he%wing(:,i), ionode_id)
          call mp_bcast(he%wing_c(:,i), ionode_id)
       enddo
    else
       nullify(he%wing)
       nullify(he%wing_c)
    endif

END SUBROUTINE read_data_pw_head_epsilon


SUBROUTINE read_data_pw_cprim_prod(cpp, prefix, l_vc, ok_read, l_vcw_overlap, l_upper)
!this subroutine read the products cprim c v\tilde{w^P} from disk
    USE kinds,                ONLY : DP
    USE io_global,            ONLY : stdout, ionode, ionode_id
    USE io_files,             ONLY : find_free_unit
    USE basic_structures,     ONLY : cprim_prod,free_memory
    USE mp,                   ONLY : mp_bcast, mp_barrier

    implicit none


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
                inquire(file=trim(prefix)//'.cprim.'//nfile,exist=ok_read)
             else
                inquire(file=trim(prefix)//'.vcprim.'//nfile,exist=ok_read)
             endif
          else
             if(.not. l_vc) then
                inquire(file=trim(prefix)//'.cprim_up.'//nfile,exist=ok_read)
             else
                inquire(file=trim(prefix)//'.vcprim_up.'//nfile,exist=ok_read)
             endif
          endif
       endif
       call mp_bcast(ok_read, ionode_id)
       if(.not. ok_read) return
    endif

    if(ionode) then
       iunsterms =  find_free_unit()
       write(nfile,'(4i1)') &
            & cpp%cprim/1000,mod(cpp%cprim,1000)/100,mod(cpp%cprim,100)/10,mod(cpp%cprim,10)
       if(.not.l_upper) then
          if(l_vcw_overlap) then
             open( unit= iunsterms, file=trim(prefix)//'.vcw_overlap.'//nfile, status='old',form='unformatted')
          else
             if(.not. l_vc) then
                open( unit= iunsterms, file=trim(prefix)//'.cprim.'//nfile, status='old',form='unformatted')
             else
                open( unit= iunsterms, file=trim(prefix)//'.vcprim.'//nfile, status='old',form='unformatted')
             endif
          endif
       else
           if(l_vcw_overlap) then
             open( unit= iunsterms, file=trim(prefix)//'.vcw_up_overlap.'//nfile, status='old',form='unformatted')
          else
             if(.not. l_vc) then
                open( unit= iunsterms, file=trim(prefix)//'.cprim_up.'//nfile, status='old',form='unformatted')
             else
                open( unit= iunsterms, file=trim(prefix)//'.vcprim_up.'//nfile, status='old',form='unformatted')
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
    call mp_bcast(cpp%nums_occ, ionode_id)
    call mp_bcast(cpp%nums, ionode_id)
    call mp_bcast(cpp%numpw, ionode_id)

    cpp%nums_cond=cpp%nums-cpp%nums_occ
    if(.not.l_vc .or. l_vcw_overlap .and. .not.l_upper) then
       allocate(cpp%cpmat(cpp%numpw,cpp%nums_cond))
    else
       allocate(cpp%cpmat(cpp%numpw,cpp%nums))
    endif
    cpp%lda=cpp%numpw
    if(.not. l_vc .or. l_vcw_overlap .and. .not.l_upper) then
       do i=1,cpp%nums_cond
          call mp_barrier
          if(ionode) read(iunsterms) cpp%cpmat(1:cpp%numpw,i)
          call mp_bcast(cpp%cpmat(:,i), ionode_id)
       enddo
    else
        do i=1,cpp%nums
           call mp_barrier
          if(ionode) read(iunsterms) cpp%cpmat(1:cpp%numpw,i)
          call mp_bcast(cpp%cpmat(:,i), ionode_id)
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
    USE io_files,             ONLY : find_free_unit
    USE basic_structures
    USE mp,                   ONLY : mp_bcast

    implicit none

    REAL(kind=DP) :: ene_dft_xc(max_i)
    INTEGER :: max_i
    CHARACTER(LEN=256) ::  prefix!to designate the PW files

    INTEGER :: iunu, i,nn


    if(ionode) then
       iunu = find_free_unit()
       open( unit=iunu, file=trim(prefix)//'.dft_xc', status='old',form='unformatted')
       read(iunu) nn
       do i=1,max_i
          read(iunu) ene_dft_xc(i)
       enddo
    endif
    call mp_bcast(ene_dft_xc(1:max_i), ionode_id)

END SUBROUTINE read_data_pw_dft_xc


SUBROUTINE read_data_pw_dft_u(ene_dft_u,max_i,prefix)
!this subroutine reads in the expectation value dur to the u energies

    USE kinds,                ONLY : DP
    USE io_global,            ONLY : stdout, ionode, ionode_id
    USE io_files,             ONLY : find_free_unit
    USE basic_structures
    USE mp,                   ONLY : mp_bcast

    implicit none

    REAL(kind=DP), INTENT(out) :: ene_dft_u(max_i)
    INTEGER, INTENT(in) :: max_i
    LOGICAL :: lex
    CHARACTER(LEN=256), INTENT(in) ::  prefix!to designate the PW files

    INTEGER :: iunu, i,nn


    if(ionode) then
       !
       INQUIRE ( file=trim(prefix)//'.hubbard_u', EXIST=lex )
       write(stdout,*) 'hubbard_u files exist ? lex=', lex
       call flush_unit(stdout)
       !
       if ( lex ) then
         !
         iunu = find_free_unit()
         open( unit=iunu, file=trim(prefix)//'.hubbard_u', status='old',form='unformatted')
         read(iunu) nn
         do i=1,max_i
            read(iunu) ene_dft_u(i)
         enddo
         !
       else
         !
         do i=1,max_i
            ene_dft_u(i)=0.0
         enddo
         !
       endif
    endif
    call mp_bcast(ene_dft_u(1:max_i), ionode_id)

END SUBROUTINE read_data_pw_dft_u


SUBROUTINE read_data_pw_upper_states(us,prefix)
!this subroutine reads in the upper REDUCED states

    USE kinds,                ONLY : DP
    USE io_global,            ONLY : stdout, ionode, ionode_id
    USE io_files,             ONLY : find_free_unit
    USE basic_structures
    USE mp,                   ONLY : mp_bcast

    implicit none

    TYPE(upper_states) :: us!structure to be read and initialized
    CHARACTER(LEN=256) ::  prefix!to designate the PW files

    INTEGER :: iunu
    INTEGER :: ii

    if(ionode) then
       iunu = find_free_unit()
       open( unit=iunu, file=trim(prefix)//'.upper', status='old',form='unformatted')
       read(iunu) us%nums_tot
       read(iunu) us%nums
       read(iunu) us%nums_occ
       read(iunu) us%nums_reduced
    endif

    call mp_bcast(us%nums_tot, ionode_id)
    call mp_bcast(us%nums, ionode_id)
    call mp_bcast(us%nums_occ, ionode_id)
    call mp_bcast(us%nums_reduced, ionode_id)

    allocate(us%ene(us%nums_reduced))

    if(ionode) then
       do ii=1,us%nums_reduced
          read(iunu) us%ene(ii)
       enddo
       close(iunu)
    endif
    call mp_bcast(us%ene(:),ionode_id)

END SUBROUTINE read_data_pw_upper_states

