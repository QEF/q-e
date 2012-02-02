! FOR GWW
! Author: P. Umari
! Modified by G. Stenuit
!
!----------------------------------------------------------------------------
SUBROUTINE produce_wannier_gamma
  !----------------------------------------------------------------------------
  !
  ! ... This routine is a driver of the self-consistent cycle.
  ! ... It uses the routine c_bands for computing the bands at fixed
  ! ... Hamiltonian, the routine sum_band to compute the charge
  ! ... density, the routine v_of_rho to compute the new potential
  ! ... and the routine mix_rho to mix input and output charge densities
  ! ... It prints on output the total energy and its decomposition in
  ! ... the separate contributions.
  !
  !
  USE kinds,                ONLY : DP
  USE constants,            ONLY : rytoev
  USE io_global,            ONLY : stdout, ionode, ionode_id
  USE wvfct,                ONLY : nbnd, et, wg,npwx,npw
  USE ener,                 ONLY : etot, eband, deband, ehart, vtxc, etxc, &
                                   etxcc, ewld, demet, ef, ef_up, ef_dw
  USE io_files,             ONLY : prefix, iunwfc, nwordwfc, find_free_unit
  USE wavefunctions_module, ONLY : evc
#if defined (EXX)
  USE funct,                ONLY : dft_is_hybrid, exx_is_active
#endif
  USE ldaU,                 ONLY : lda_plus_u
  USE wannier_gw,           ONLY : u_trans,numw_prod, real_matrix_pointer, complex_matrix_pointer, &
                                 &    wannier_P, w_prod, free_memory, num_nbndv, &
                                 &    cutoff_overlap, nset, num_nbnds, num_nbndc_set, expgsave, &
                                 &    remainder, restart_gww, num_nbnd_first, l_zero, l_wing, &
                                 &    cprim_type, l_vcw_overlap, nset_overlap, nspace, lambda_ene, l_orthonorm_products, &
                                 &     cutoff_products, ecutoff_global,numw_prod_vvc,&
                                 &    l_wpwt_terms, l_polarization_analysis, &
                                 &    cutoff_polarization, nc_polarization_analysis, l_only_val_cond, numw_prod_val_cond, &
                                 &    l_no_val_cond_sec, numw_prod_val_cond_sec, l_plot_mlwf, l_plot_ulwf, l_ultra_external,&
                                 &    no_radius, nbnd_normal, max_ngm_set, l_pmatrix, p_mpime, p_nproc, npcol, nprow, icontxt, &
                                 &    myrow, mycol, l_coulomb_analysis,cutoff_coulomb_analysis, truncation_radius
  USE uspp,                 ONLY : okvan
  USE mp,                   ONLY : mp_bcast,mp_barrier
  USE mp_global,            ONLY : nproc, mpime
  USE parallel_include
  USE exx,                  ONLY : exx_grid_init
  USE buffers,              ONLY : get_buffer, save_buffer
  !
  IMPLICIT NONE
  !
  ! ... a few local variables
  !
#if defined (EXX)
  REAL(DP) :: dexx
  REAL(DP) :: fock0,  fock1,  fock2
#endif
  REAL(kind=DP), TARGET, ALLOCATABLE :: omat(:,:)!matrix for overlap between product of wanniers
  REAL(kind=DP), TARGET, ALLOCATABLE :: uterms(:,:), uterms_tmp(:)!matrix for 1/|r-r'| terms between product of wanniers
  TYPE(wannier_P), ALLOCATABLE :: w_P(:)
  INTEGER :: iunuterms, iw, jw
  INTEGER :: numw_prod_old
  REAL(kind=DP), ALLOCATABLE :: vmat(:,:)
  INTEGER :: numpw_max, numpw_offset
  !
  REAL(kind=DP), ALLOCATABLE :: e_xc(:), e_h(:), e_u(:)
  INTEGER :: ii
  REAL(kind=DP), ALLOCATABLE :: ene_loc(:)!local energy for state
  LOGICAL :: l_on_ene!if true localizes also mon energy
  REAL(kind=DP), ALLOCATABLE :: tmp_rot(:,:)
  INTEGER :: numw_prod_r, numw_prod_c, numw_prod_dimr, numw_prod_dimc
  INTEGER :: num_wp_r, num_wp_dimr
#ifdef __SCALAPACK
  INTEGER, EXTERNAL :: indxg2p, indxg2l
#endif
  INTEGER :: numw_prod_vvc_r, numw_prod_vvc_dimr, numw_prod_vvc_c, numw_prod_vvc_dimc

  INTEGER :: i

! NEW
  allocate(e_xc(nbnd),e_h(nbnd))
  if (lda_plus_u) allocate(e_u(nbnd))

!setup parallel environment
#ifndef __MPI
         l_pmatrix=.false.
#endif
#ifndef __SCALAPACK
         l_pmatrix=.false.
#endif

         if(l_pmatrix) then
#ifdef __SCALAPACK
            call blacs_pinfo(p_mpime,p_nproc)
            write(stdout,*) 'PINFO',p_mpime,p_nproc
       !     nprow=int(sqrt(real(p_nproc)))
       !     npcol=p_nproc/nprow
            write(stdout,*) 'NPROW NPCOL', nprow, npcol

            call blacs_get(0,0,icontxt)
            call blacs_gridinit(icontxt,'R',nprow, npcol)
            call blacs_gridinfo(icontxt, nprow, npcol, myrow,mycol)
            write(stdout,*) 'MYROW MYCOL', myrow,mycol
#endif
         endif



         !set ngm max
         call max_ngm_set

         if(lambda_ene > 0 .and. lambda_ene <= 1.d0) then
            l_on_ene = .true.
         else
            l_on_ene = .false.
         endif

         allocate(ene_loc(nbnd))
! ene_loc is defined as inout in go_wannier and is re-used for restart_gw<=3 in go_wannier. NO
! cprim_type allways=2 (default value)
         ene_loc(:)= et(:,1)*rytoev
         call exx_grid_init()

!         if(remainder /= 4 .and. remainder /= 5 )then !otherwise just post-processing remainder calculation
            !
            if(restart_gww <= 0) then
               write(stdout,*) 'restart_gww <= 0'
               write(stdout,*) "remainder", remainder
               write(stdout,*) 'ATT1'
               write(stdout,*) 'evc(1,1)=', evc(1,1), 'evc(1,2)=', evc(1,2)
               write(stdout,*) npwx, npw, nbnd, evc(1,1)
               write(stdout,*) rytoev
               do i=1,nbnd
                 write(stdout,*)i,et(i,1)*rytoev
               enddo
               call flush_unit(stdout)
               CALL wfc_gamma_real(0)
               write(stdout,*) 'ATT2'
               call flush_unit(stdout)
               do i=1,nbnd
                 write(stdout,*)i,et(i,1)*rytoev
               enddo
               write(stdout,*) npwx, npw, nbnd, evc(1,1)
               !!!! only 5 arguments !!!!
               !!! compared to the previous version, the 4th
               !!! argument, evc, has been removed !!!
               call energies_xc( npwx, npw, nbnd, e_xc, e_h )
               write(stdout,*) 'After : call energies_xc'
               write(stdout,*) '------------------------------'
               write(stdout,*) 'e_xc(1:nbnd)=', e_xc(1:nbnd)
               write(stdout,*) '------------------------------'
               write(stdout,*) 'e_h(1:nbnd)=', e_h(1:nbnd)
               write(stdout,*) '------------------------------'
               CALL flush_unit( stdout )
               !
               call write_energies_xc(e_xc)
               !
               write(stdout,*) 'After : call write_energies_xc'
               write(stdout,*) '------------------------------'
               write(stdout,*) 'e_xc(1:nbnd)=', e_xc(1:nbnd)
               write(stdout,*) '------------------------------'
               write(stdout,*) 'e_h(1:nbnd)=', e_h(1:nbnd)
               write(stdout,*) '------------------------------'
               CALL flush_unit( stdout )

               !!! NEW:
               IF ( lda_plus_u ) then
                  write(stdout,*) 'if lda_plus_u : compute <ener>_i'
                  !!!write(stdout,*) npwx, npw, nbnd, evc(1,1)
                  CALL flush_unit( stdout )
                  CALL energies_u_gamma( e_u)
                  write(stdout,*) 'write to a file .hubbard_u these energies'
                  !!!write(stdout,*) npwx, npw, nbnd, evc(1,1)
                  DO iw=1, nbnd
                     write(stdout,*) 'e_u(', iw , ')=', e_u(iw), ' Ry'
                  ENDDO
                  CALL flush_unit( stdout )
                  CALL write_energies_u(e_u)
                  write(stdout,*) 'OUT of lda_plus_u'
                  CALL flush_unit( stdout )
               endif

               write(stdout,*) 'call go_wannier'
               CALL flush_unit( stdout )
               !-----------------------------------------------
               !----------------------------------------------
               call go_wannier(iunwfc,1.d-9,40,num_nbndv, 0, l_on_ene, ene_loc, lambda_ene)
               !
               !!!! just before, in go_wannier, the evc have been changed
               !!!! by rotate_wannier_gamma and stored in evc1 in files wfc_w
               !!!! so the evc have not been updated !!!! => one need here to read back the evc
               !!!! or changed wfc_gamma_real !!!
!!!!!               call wfc_gamma_real(0)
! MODIFIED to
               call wfc_gamma_real_after_rot(0)
               ! reput it to check if QE will be OK
! pourquoi????               call wfc_gamma_real(0)
               !if required save MLWF for plotting
               write(stdout,*) 'call write_wfc_plot if l_plot_mlwf .or. l_ultra_external is T'
               CALL flush_unit( stdout )
               if(l_plot_mlwf.or.l_ultra_external) then
                  call write_wfc_plot(0)
               endif
               !               deallocate(evc)
               call distance_wannier
               CALL flush_unit( stdout )
               ! call appropriate routine for ultralocalization
               if(.not.l_ultra_external) then
                  call ultralocalization_para(num_nbndv,nbnd_normal,1.d-5,0,num_nbndv,0)
                  write(stdout,*) 'Ultralocalization valence wfcs'!ATTENZIONE
                  if(num_nbnd_first==0) then
                     call ultralocalization_para(num_nbndv,nbnd_normal,1.d-5,1,nbnd-num_nbndv,0)
                  else
                     call ultralocalization_para(num_nbndv,num_nbndv+num_nbnd_first,1.d-5,1,num_nbnd_first,0)
                     call ultralocalization_para(num_nbndv+num_nbnd_first,nbnd_normal,1.d-5,2,&
                          &nbnd-num_nbndv- num_nbnd_first,0)
                  endif
               else
                   call ultra_external( 1, num_nbndv, no_radius, 0)
                   if(num_nbnd_first==0) then
                      call ultra_external( num_nbndv+1, nbnd_normal, no_radius, 0)
                   else
                      call ultra_external( num_nbndv+1, num_nbndv+num_nbnd_first, no_radius, 0)
                      call ultra_external( num_nbndv+num_nbnd_first+1, nbnd_normal, no_radius, 0)
                   endif
               endif
               if(okvan) deallocate(expgsave)
               !if required save MLWF for plotting
               if(l_plot_ulwf) then
                  !                   ALLOCATE( evc( npwx, nbnd ) )
                  !                   CALL get_buffer ( evc, nwordwfc, iunwfc, 1)
                  allocate(tmp_rot(nbnd_normal,nbnd_normal))
                  tmp_rot(:,:)=dble(u_trans(:,:))
                  call rotate_wannier_gamma( tmp_rot,1,1)
                  deallocate(tmp_rot)
                  call write_wfc_plot(1)
                  !                 deallocate(evc)
               endif
               !calculate products of wanniers


               !write transformation matrix u on file
               !!!
               ! to check
               write(stdout,*) '------------------------------'
               write(stdout,*) 'e_xc(1:nbnd)=', e_xc(1:nbnd)
               write(stdout,*) '------------------------------'
               write(stdout,*) 'e_h(1:nbnd)=', e_h(1:nbnd)
               write(stdout,*) '------------------------------'
               CALL flush_unit( stdout )
               !
               if(ionode) call write_wannier_matrix(e_xc,e_h)
               do ii=1,nbnd
                  call mp_barrier
                  call mp_bcast(u_trans(:,ii),ionode_id)
               enddo
               ! u_trans TO BE BROADCASTED
               deallocate(e_xc,e_h)
               !
               call product_wannier_para(num_nbndv,.true., ene_loc, lambda_ene)
               !
               ! now it can delete the wavefunctions in real space
               if(cprim_type /= 0) then
                  call dirdel('real_whole')
                  call dirdel('realwan')
               endif
               !
               CALL flush_unit( stdout )
               ! calculate terms for calculating residual part
               ! not in use now
               !
               !               if(remainder==0) then
               !                  call wannier_valence_terms(num_nbnds,nset,num_nbnds)
               !               else if(remainder==1) then
               !                  call wannier_valence_terms_cutoff(num_nbnds,nset,num_nbnds)
               !               else if(remainder==3) then
               !                  call wannier_valence_terms_distance(num_nbnds,nset,num_nbnds)
               !               endif
               !               deallocate(l_on_products)!ATTENZIONE posizione provvisoria
               !
               !calculate 1/|r-r'| between products of wannier
            else
               !                deallocate(evc)
            endif  !!!!!!! endif gww<=0   !!!!! endif of the first restart gww
            !
            if(restart_gww <= 1 ) then
                write(stdout,*) 'restart_gww <= 1 and numw_prod= ', numw_prod
                call flush_unit(stdout)
                !
                if(l_orthonorm_products) then
                   if(.not.l_pmatrix) then
                      allocate(omat(numw_prod,numw_prod))
                      numw_prod_dimr=numw_prod
                      numw_prod_dimc=numw_prod
                      numw_prod_r=numw_prod
                      numw_prod_c=numw_prod
                   else
#ifdef __SCALAPACK
                      numw_prod_r=ceiling(real(numw_prod)/real(max(nprow,npcol)))
                      numw_prod_c=ceiling(real(numw_prod)/real(max(npcol,nprow)))
                      numw_prod_dimr=ceiling (real(numw_prod)/real(numw_prod_r*nprow))*numw_prod_r
                      numw_prod_dimc=ceiling (real(numw_prod)/real(numw_prod_c*npcol))*numw_prod_c
                      allocate(omat(numw_prod_dimr,numw_prod_dimc))

#endif
                   endif
                   write(stdout,*) 'CALL wannier_pmat_terms_ggrid',numw_prod_dimr,numw_prod_dimc
                   call flush_unit(stdout)
                   if(.not.l_pmatrix) then
                      call wannier_pmat_terms_ggrid(numw_prod,numw_prod,omat, nset, ecutoff_global)
                   else
                      call wannier_pmat_terms_ggrid(numw_prod_dimr,numw_prod_dimc,omat, nset, ecutoff_global)
                   endif
                   write(stdout,*) 'CALL orthonormalize_producs_cutoff', numw_prod_dimr,numw_prod_dimc
                   write(stdout,*) 'OMAT1', omat(1,1)
                   call flush_unit(stdout)
                   numw_prod_vvc=numw_prod
                   !
                   if(.not.l_only_val_cond) then
                      numpw_max=numw_prod
                   else
                      numpw_max=numw_prod_val_cond
                   endif
                   !
                   if(.not.l_no_val_cond_sec) then
                      numpw_offset=0
                   else
                      numpw_offset=numw_prod_val_cond_sec
                   endif
                   !
                   if(restart_gww==1) then
                   write(stdout,*) 'restart_gww == 1 and numw_prod= ', numw_prod
                      call flush_unit(stdout)
                      !we suppose ordinary alpha startegy calculation
                      numpw_offset=0
                      numpw_max=numw_prod
                   endif
                   !
                   if(.not.l_pmatrix) then
                      call  orthonormalize_producs_cutoff(numw_prod,numpw_max, &
                      numpw_offset,numw_prod_vvc,omat,numw_prod,numw_prod,cutoff_products, nset, &
                           & ecutoff_global,.true.)
                   else
                      call  orthonormalize_producs_cutoff(numw_prod,numpw_max,&
                      numpw_offset,numw_prod_vvc,omat,numw_prod_dimr,numw_prod_dimc,cutoff_products, nset,&
                           & ecutoff_global,.true.)
                   endif
                   !
                   write(stdout,*) 'OMAT2', omat(1,1)
                   !
                   if(l_wpwt_terms.or.l_polarization_analysis) then
                     !
                     !find the overlaps <\tilde{w}_i*\tilde{}w_j|w^P_k>
                     if(.not.l_pmatrix) then
                        allocate(w_P(numw_prod))
                        num_wp_r=numw_prod
                        num_wp_dimr=numw_prod
                     else
                        num_wp_r=ceiling(real(numw_prod)/real(nproc))
                        num_wp_dimr=ceiling (real(numw_prod)/real(numw_prod_r*nproc))*num_wp_r
                        allocate(w_P(num_wp_dimr))
                     endif
                     write(stdout,*) 'CALL set_wannier_P'
                     call flush_unit(stdout)
                     !
                     !
                     call set_wannier_P(numw_prod,numw_prod_vvc,w_prod,omat,numw_prod_dimr,&
                     numw_prod_dimc,cutoff_overlap,w_P, num_wp_dimr, num_wp_r, numw_prod_r,numw_prod_c)
                     CALL flush_unit( stdout )
                     !
                     ! write overlaps  <w_P|w_i*w_j> on file
                     call write_wannier_products(w_P,0,num_wp_dimr,num_wp_r)
                     !
                     write(stdout,*) 'OMAT3', omat(1,1)
                     !write orthonormalization matrix of wanniers products
                     !
                     numw_prod_old=numw_prod
                     ! if required do the polarizarion analysis
                     if(l_polarization_analysis) then
                        call do_polarization_analysis( w_P,num_wp_dimr,num_wp_r, &
                        cutoff_polarization,numw_prod,&
                              &numw_prod_vvc, nset, ecutoff_global, .true.,omat,numw_prod_dimr,&
                              numw_prod_dimc,nc_polarization_analysis,numw_prod_r,numw_prod_c)
                        !
                        write(stdout,*) 'OMAT4', omat(1,1)
                        !
                        ! deallocate memory
                        !ATTENZIONE the following is now done in do_polarization_analysis
                        !                         do iw=1,numw_prod_old
                        !                            call free_memory(w_P(iw))
                        !                         enddo
                        deallocate(w_P)
                        !
                        if(l_coulomb_analysis) then
                           !allocate and create coulomb matrix
                           if(.not.l_pmatrix) then
                              numw_prod_r=numw_prod
                              numw_prod_dimr=numw_prod
                              numw_prod_c=numw_prod
                              numw_prod_dimc=numw_prod
                           else
                              numw_prod_r=ceiling(real(numw_prod)/real(max(nprow,npcol)))
                              numw_prod_c=ceiling(real(numw_prod)/real(max(npcol,nprow)))
                              numw_prod_dimr=ceiling (real(numw_prod)/real(numw_prod_r*nprow))*numw_prod_r
                              numw_prod_dimc=ceiling (real(numw_prod)/real(numw_prod_c*npcol))*numw_prod_c
                           endif
                           !
                           if(.not.l_pmatrix) then
                              numw_prod_vvc_r=numw_prod_vvc
                              numw_prod_vvc_dimr=numw_prod_vvc
                              numw_prod_vvc_c=numw_prod_vvc
                              numw_prod_vvc_dimc=numw_prod_vvc
                           else
                              numw_prod_vvc_r=ceiling(real(numw_prod_vvc)/real(max(nprow,npcol)))
                              numw_prod_vvc_c=ceiling(real(numw_prod_vvc)/real(max(npcol,nprow)))
                              numw_prod_vvc_dimr=ceiling (real(numw_prod_vvc)/ &
                              real(numw_prod_vvc_r*nprow))*numw_prod_vvc_r
                              numw_prod_vvc_dimc=ceiling (real(numw_prod_vvc)/&
                              real(numw_prod_vvc_c*npcol))*numw_prod_vvc_c
                           endif
                           !
                           allocate(vmat(numw_prod_dimr,numw_prod_dimc))
                           !calculate vmat matrix
                           call  wannier_uterms_red(nset,l_zero, ecutoff_global,vmat,numw_prod_dimr,numw_prod_dimc,&
                                 &numw_prod_r,numw_prod_c,numw_prod,numw_prod,1,.true.)

                           !call routine
                           call coulomb_analysis(numw_prod,numw_prod_old,numw_prod_vvc,vmat, numw_prod_dimr, &
                                  &numw_prod_dimc,cutoff_coulomb_analysis, nset, ecutoff_global, .true.,omat, &
                                  &numw_prod_vvc_dimr,numw_prod_vvc_dimc,numw_prod_vvc_r,numw_prod_vvc_c)
                           !
                           !deallocate
                           deallocate(vmat)
                        endif  !!! of if(l_coulomb_analysis)
                        ! set new descriptors for products of wanniers
                        if(.not.l_pmatrix) then
                           allocate(w_P(numw_prod))
                           num_wp_r=numw_prod
                           num_wp_dimr=numw_prod
                        else
                           num_wp_r=ceiling(real(numw_prod)/real(nproc))
                           num_wp_dimr=ceiling (real(numw_prod)/real(numw_prod_r*nproc))*num_wp_r
                           allocate(w_P(num_wp_dimr))
                        endif
                        !
                        call set_wannier_P(numw_prod,numw_prod_vvc,w_prod,omat,numw_prod_dimr,&
                        numw_prod_dimc,cutoff_overlap,w_P,num_wp_dimr,num_wp_r,numw_prod_r,numw_prod_c)
                        CALL flush_unit( stdout )
                        !
                        ! write overlaps  <w_P|w_i*w_j> on file
                        !
                        call write_wannier_products(w_P,0,num_wp_dimr,num_wp_r)
                        if(.not.l_pmatrix) then
                           do iw=1,numw_prod
                              if(ionode) call free_memory(w_P(iw))
                           enddo
                        else
#ifdef __SCALAPACK
                           do iw=1,numw_prod
                              if(indxg2p(iw,num_wp_r,0,0,nproc)==mpime) then
                                 call free_memory(w_P(indxg2l(iw,num_wp_r,0,0,nproc)))
                              endif
                           enddo
#endif
                        endif
                        deallocate(w_P)
                        deallocate(omat)
                        !the following is for the beta strategy
                        if(cprim_type /= 2) then
                           !TD allocate vmat
                           if(.not.l_pmatrix) then
                              numw_prod_vvc_r=numw_prod_vvc
                              numw_prod_vvc_dimr=numw_prod_vvc
                              numw_prod_r=numw_prod
                              numw_prod_dimr=numw_prod
                              numw_prod_c=numw_prod
                              numw_prod_dimc=numw_prod
                           else
                              numw_prod_r=ceiling(real(numw_prod)/real(max(nprow,npcol)))
                              numw_prod_c=ceiling(real(numw_prod)/real(max(npcol,nprow)))
                              numw_prod_dimr=ceiling (real(numw_prod)/real(numw_prod_r*nprow))*numw_prod_r
                              numw_prod_dimc=ceiling (real(numw_prod)/real(numw_prod_c*npcol))*numw_prod_c
                              numw_prod_vvc_r=ceiling(real(numw_prod_vvc)/real(nprow))
                              numw_prod_vvc_dimr=ceiling (real(numw_prod_vvc)/real(numw_prod_vvc_r*nprow))*numw_prod_vvc_r
                           endif
                           !
                           allocate(vmat(numw_prod_vvc_dimr,numw_prod_dimc))
                           !calculate vmat matrix
                           !
                           call  wannier_uterms_red(nset,l_zero, ecutoff_global,vmat,numw_prod_vvc_dimr,numw_prod_dimc,&
                                 &numw_prod_vvc_r,numw_prod_c,numw_prod,numw_prod_vvc,0,.false.)
                           !treat is as products of wanniers
                           !
                           if(.not.l_pmatrix) then
                              allocate(w_P(numw_prod))
                              num_wp_r=numw_prod
                              num_wp_dimr=numw_prod
                           else
                              num_wp_r=ceiling(real(numw_prod)/real(nproc))
                              num_wp_dimr=ceiling (real(numw_prod)/real(num_wp_r*nproc))*num_wp_r
                              allocate(w_P(num_wp_dimr))
                           endif
                           call set_wannier_P(numw_prod,numw_prod_vvc,w_prod,vmat,numw_prod_vvc_dimr,numw_prod_dimc,&
                                 &0.d0,w_P,num_wp_dimr,num_wp_r,numw_prod_vvc_r,numw_prod_c)

                           !write on disk
                           call write_wannier_products(w_P,1,num_wp_dimr,num_wp_r)
                           !
                           if(.not.l_pmatrix) then
                               do iw=1,numw_prod
                                  if(ionode)call free_memory(w_P(iw))
                               enddo
                           else
#ifdef __SCALAPACK
                               do iw=1,numw_prod
                                  if(indxg2p(iw,num_wp_r,0,0,nproc)==mpime) then
                                     call free_memory(w_P(indxg2l(iw,num_wp_r,0,0,nproc)))
                                  endif
                               enddo
#endif
                           endif
                           deallocate(w_P)
                           deallocate(vmat)
                        else !!! if(cprim_type /= 2)
                            !call dirdel('wiwjwfc')!ATTENZIONE
                        endif
                        call dirdel('wiwjwfc_on')
                     else  !!! f(l_polarization_analysis)
                        if(.not. l_pmatrix) then
                           do iw=1,numw_prod_old
                              if(ionode) call free_memory(w_P(iw))
                           enddo
                        else
#ifdef __SCALAPACK
                          do iw=1,numw_prod_old
                             if(indxg2p(iw,num_wp_r,0,0,nproc)==mpime) then
                                 call free_memory(w_P(indxg2l(iw,num_wp_r,0,0,nproc)))
                             endif
                          enddo
#endif
                        endif !!!! endif of if(.not. l_pmatrix)
                        deallocate(omat)
                     endif !!!!!!! endif if(l_polarization_analysis)
                  else  !!!!! else of if(l_wpwt_terms.or.l_pola ...
                     deallocate(omat)
                  endif !!!!! endif of if if(l_wpwt_terms.or.l_pola ...
                endif  !!!! if(l_orthonorm_products)
                !
                if(.not. l_zero) then
                   if(.not.l_orthonorm_products) then
                      call wannier_uterms(nset,.false.,.false.,0, ecutoff_global)
                   else
                      if(.not.l_polarization_analysis) then
                         call wannier_uterms(nset,.false.,.false.,1, ecutoff_global)
                      else
                         call wannier_uterms(nset,.false.,.false.,2, ecutoff_global)
                      endif
                   endif
                   CALL flush_unit( stdout )
                   allocate(uterms(numw_prod,numw_prod))
                   if(ionode) then
                      iunuterms =  find_free_unit()
                      open( unit= iunuterms, file=trim(prefix)//'.uterms', status='old',form='unformatted')
                   endif
                   allocate(uterms_tmp(numw_prod))
                   do iw=1,numw_prod
                      if(ionode) read(iunuterms) uterms_tmp(1:iw)
                      call mp_bcast(uterms_tmp(1:iw), ionode_id)
                      do jw=1,iw
                         uterms(iw,jw)=uterms_tmp(jw)
                         uterms(jw,iw)=uterms(iw,jw)
                      enddo
                   enddo
                   deallocate(uterms_tmp)
                   if(ionode)    close(iunuterms)
                   if(ionode) call write_vpot_matrix(uterms,0)
                   deallocate(uterms)
                   CALL flush_unit( stdout )
                   !
                   !if required also the term with G=0,G'=0 of v put to zero
                else !!!! if(.not. l_zero)
                   if(.not.l_orthonorm_products) then
                       call wannier_uterms(nset,.false.,.true.,0, ecutoff_global)
                   else
                      if(.not.l_polarization_analysis) then
                         call wannier_uterms(nset,.false.,.true.,1, ecutoff_global)
                      else
                         call wannier_uterms(nset,.false.,.true.,2, ecutoff_global)
                      endif
                   endif
                   CALL flush_unit( stdout )
                   allocate(uterms(numw_prod,numw_prod))
                   if(ionode) then
                     iunuterms =  find_free_unit()
                     open( unit= iunuterms, file=trim(prefix)//'.uterms', status='old',form='unformatted')
                   endif
                   write(stdout,*) 'WRITE UTERMS'
                   call flush_unit(stdout)
                   call mp_barrier
                   !
                   allocate(uterms_tmp(numw_prod))
                   do iw=1,numw_prod
                     if(ionode) read(iunuterms) uterms_tmp(1:iw)
                     call mp_barrier
                     call mp_bcast(uterms_tmp(1:iw), ionode_id)
                     do jw=1,iw
                       uterms(iw,jw)=uterms_tmp(jw)
                       uterms(jw,iw)=uterms(iw,jw)
                     enddo
                   enddo
                   deallocate(uterms_tmp)
                   if(ionode)    close(iunuterms)
                   if(ionode) call write_vpot_matrix(uterms,3)
                   deallocate(uterms)
                   CALL flush_unit( stdout )
                endif !!!!! !!!! endif of if(.not. l_zero)
                !
                write(stdout,*) 'OUT OF RESTART_GWW1',numw_prod
                call flush_unit(stdout)
                call mp_barrier
                !
            endif  !!!! endif restartgw<=1

            if(restart_gww <= 2) then
              !!!! for after ...
            endif
            !
            if(restart_gww <= 3 ) then
              write(stdout,*) 'restart_gww <= 3 and numw_prod= ', numw_prod
              call flush_unit(stdout)
              !
              !if required calculates the wings term of the symmetrized dielectric matrix
              if( l_wing) then
                 if(.not.l_orthonorm_products) then
                    call calculate_wing(nset,0)
                 else
                    if(.not.l_polarization_analysis) then
                       call calculate_wing(nset,1)
                    else
                       call calculate_wing(nset,2)
                    endif
                 endif
              endif
              !
              ! if required ultralocalized subspace {C'} and calculates Vc',c'p matrix elements
              if(cprim_type == 0) then !usual treatment
                 if(num_nbndc_set > 0 ) then
                    !
                    !                     ALLOCATE( evc( npwx, nbnd ) )
                    !                     CALL get_buffer ( evc, nwordwfc, iunwfc, 1)
                    ! pourquoi???         CALL wfc_gamma_real(1)

                    !----------------------------------------------
                    ! I put here ene_loc to re-initialize
!                    ene_loc(:)= et(:,1)*rytoev
                    !-----------------------------------------------
                    call go_wannier(iunwfc,1.d-9,40,num_nbndv, 1,.false., ene_loc, 0.1d0)
                    ! pourquoi???                     call wfc_gamma_real(1)
                    !if required save MLWF for plotting
                    if(l_plot_mlwf.or.l_ultra_external) then
                        call write_wfc_plot(0)
                    endif
!                   deallocate(evc)
                    if(.not.l_ultra_external) then
                       call ultralocalization_para(num_nbndv,nbnd_normal,1.d-5,1,nbnd-num_nbndv,1)
                    else
                       call ultra_external( num_nbndv+1, num_nbndv+num_nbndc_set, no_radius, 1)
                    endif
                    if(restart_gww == 0) then
                        call product_wannier_para_c(num_nbndv,.true.,.false.)
                     else
                        call product_wannier_para_c(num_nbndv,.true.,.true.)
                     endif
                     if(.not.l_zero) then
                        if(.not.l_orthonorm_products) then
                           call wannier_uterms_c(nset,.false.,0, ecutoff_global)
                        else
                           if(.not.l_polarization_analysis) then
                              call wannier_uterms_c(nset,.false.,1, ecutoff_global)
                           else
                              call wannier_uterms_c(nset,.false.,2, ecutoff_global)
                           endif
                        endif
                        CALL flush_unit( stdout )
!if required calculated also the term with v(G=0,G'=0) = 0
                     else
                        if(.not.l_orthonorm_products) then
                           call wannier_uterms_c(nset,.true.,0, ecutoff_global)
                        else
                           if(.not.l_polarization_analysis) then
                              call wannier_uterms_c(nset,.true.,1, ecutoff_global)
                           else
                              call wannier_uterms_c(nset,.true.,2, ecutoff_global)
                           endif
                        endif
                        CALL flush_unit( stdout )

                     endif
                     call write_wannier_matrix_c



                  endif
              else if (cprim_type==1) then!direct treatment of selected states
!read in bands again
!                    ALLOCATE( evc( npwx, nbnd ) )
!                    CALL get_buffer ( evc, nwordwfc, iunwfc, 1)
                    if(.not.l_orthonorm_products) then
                       call create_vcprim(nset_overlap, l_zero,0, ecutoff_global,.true.)
                    else
                       if(.not.l_polarization_analysis) then
                          call create_vcprim(nset_overlap, l_zero,1, ecutoff_global,.true.)
                       else
                          call create_vcprim(nset_overlap, l_zero,2, ecutoff_global,.true.)
                       endif
                    endif

!                    DEALLOCATE( evc)

                 else
!direct treatment of selected states in bands again
                    write(stdout,*) 'Inside'
                    write(stdout,*) 'evc(1,1)=', evc(1,1), 'evc(1,2)=', evc(1,2)
                    call flush_unit(stdout)
                    !!!! since it seems the evc have changed ! I have removed these comments
!                    ALLOCATE( evc( npwx, nbnd ) )
                    !evc(:,:) = (0.D0,0.D0)
!                    CALL davcio(evc, nwordwfc, iunwfc, 1, - 1)
!                   ALLOCATE( evc( npwx, nbnd ) )
!                   CALL get_buffer ( evc, nwordwfc, iunwfc, 1)
                    write(stdout,*) 'after get_buffer'
                    write(stdout,*) 'evc(1,1)=', evc(1,1), 'evc(1,2)=', evc(1,2)
                    call flush_unit(stdout)
                    if(.not.l_orthonorm_products) then
                       call create_vcprim(nset_overlap, l_zero,0, ecutoff_global,.false.)
                    else
                       if(.not.l_polarization_analysis) then
                           call create_vcprim(nset_overlap, l_zero,1, ecutoff_global,.false.)
                        else
                           call create_vcprim(nset_overlap, l_zero,2, ecutoff_global,.false.)
                        endif
                     endif
!                   DEALLOCATE( evc)
                 endif
            endif
!            if(restart_gww <= 4  ) then
!               write(stdout,*) 'restart_gww <= 4 and numw_prod= ', numw_prod
!               call flush_unit(stdout)

!               CALL flush_unit( stdout )
!!               ALLOCATE( evc( npwx, nbnd ) )
!-------------------------------------------------------------
! here is just because in rotate_wannier evc have been modified.....!!!! YEAH!
!               evc(:,:) = (0.D0,0.D0)
!               CALL davcio(evc, nwordwfc, iunwfc, 1, - 1)
!!               CALL get_buffer ( evc, nwordwfc, iunwfc, 1)
!--------------------------------------------------------------
!               CALL dft_exchange(num_nbndv,num_nbnds,nset)
!               CALL flush_unit( stdout )
!!               deallocate(evc)
!            endif
            if(restart_gww <= 5  ) then
              write(stdout,*) 'restart_gww <= 5 and numw_prod= ', numw_prod
              call flush_unit(stdout)

!if required calculates the overlaps v c w_p
              CALL flush_unit( stdout )
!                 ALLOCATE( evc( npwx, nbnd ) )
!                 CALL get_buffer ( evc, nwordwfc, iunwfc, 1)
              if(l_vcw_overlap) then
                 if(.not.l_orthonorm_products) then
                    call create_vcw_overlap(nset_overlap,0,ecutoff_global)
                 else
                    if(.not.l_polarization_analysis) then
                        call create_vcw_overlap(nset_overlap,1,ecutoff_global)
                    else
                        call create_vcw_overlap(nset_overlap,2,ecutoff_global)
                    endif
                 endif
              endif
!                 deallocate(evc)
            endif ! endif restart_gww<=5
!if required calculates the reduced conduction states
              if(restart_gww <= 6  ) then
                 write(stdout,*) 'restart_gww <= 6 and numw_prod= ', numw_prod
                 call flush_unit(stdout)

                 CALL flush_unit( stdout )
!                 ALLOCATE( evc( npwx, nbnd ) )
!                 CALL get_buffer ( evc, nwordwfc, iunwfc, 1)
                 if(nbnd > nbnd_normal) then
                    if(.not.l_orthonorm_products) then
                       call create_upper_states(nset_overlap, l_zero, 0,ecutoff_global)
                    else
                       if(.not.l_polarization_analysis) then
                          call create_upper_states(nset_overlap, l_zero, 1,ecutoff_global)
                       else
                          call create_upper_states(nset_overlap, l_zero, 2,ecutoff_global)
                       endif
                    endif
                 endif
              endif
!          else!just study of exchange convergence
!             write(stdout,*) 'EXCHANGE RADIUS',truncation_radius
!             call flush_unit(stdout)
!             CALL dft_exchange(num_nbndv,num_nbnds,nset)
!             CALL flush_unit( stdout )
!          endif ! endif remainder
          deallocate(ene_loc)


!switch off parallel environment for matrices
          if(l_pmatrix) then
#ifdef __SCALAPACK
            if((myrow+1) <= nprow .and. (mycol+1) <= npcol) then
               call blacs_gridexit(icontxt)
            endif

#endif
         endif
         return
END SUBROUTINE  produce_wannier_gamma
