!
! Copyright (C) 2001-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!



MODULE basic_structures
!this module describes the basis structures 
!which are obtained from a DFT code
  USE kinds, ONLY : DP
  
  TYPE wannier_u
!this structure describes the  transformation
!from KS to Wannier states together with the KS eigenenergies
    INTEGER :: nspin!spin multiplicity
    INTEGER :: nums!number of states
    INTEGER :: nums_occ(2)!number of occupied states for the two spin channnels
    REAL(kind=DP), DIMENSION(:,:), POINTER :: ene!KS energies
    REAL(kind=DP), DIMENSION(:,:), POINTER :: ene_xc!LDA exchange and correlation terms
    REAL(kind=DP), DIMENSION(:,:), POINTER :: ene_lda_h!LDA exchange and correlation terms
    COMPLEX(kind=DP), DIMENSION(:,:,:), POINTER :: umat!INVERSE transformation matrix to wannier: Psi_i=U_{i,j}w_j
  END TYPE wannier_u


  TYPE wannier_P
!this structure described the localized and normalized products of wanniers  w_P
    INTEGER :: numij!number of (unordered)couples of wannier w_i w_j which overlap both with w_P
    INTEGER, DIMENSION(:,:), POINTER :: ij!array for i,j (ij,:)
    REAL(kind=DP),DIMENSION(:), POINTER :: o!overlap <w_P|w_i*w_j>
  END TYPE wannier_P

  TYPE v_pot
!this structure describes the coulomb potential on the base of (orthonormalized)
!products of wanniers
    INTEGER :: numpw!number of states(products)
    REAL(kind=DP),DIMENSION(:,:), POINTER :: vmat!potentail 1/|r-r'|
  END TYPE v_pot

 TYPE q_mat
!this structures describes the set of overlap of othonormalized products
!of wanniers with products of wannier
    INTEGER :: numpw!number of states(orthonormalized products)
!parameters used for  parallelization
    LOGICAL :: is_parallel!if true is a parallel part of the global q_mat matrix
    INTEGER :: numpw_para!numer of states(orthonormalized products) on this processor
    INTEGER :: first_para!first state (orthonormalized products, global order)on this processor
    TYPE(wannier_P), DIMENSION(:), POINTER :: wp!arrays of wannier products descriptors
 END TYPE q_mat


 TYPE ortho_polaw
!this structure describe the orthonormalization matrix
!w^P_i=A_{i,j}\tilde{w^P}_j
!it is put here because it's read from a PW file
    INTEGER :: numpw!number of states (products of wanniers)
    LOGICAL :: inverse!if true, the inverse transform is stored
    REAL(kind=DP), DIMENSION(:,:), POINTER :: on_mat!the transformation
 END TYPE ortho_polaw

 TYPE wp_psi
!this structure describe the product of KS wavefunctions with unorthonormalized
!products of wannier \int dr w^P_i(r) w^P_j(r) Psi_v(r) Psi_v(r)
    INTEGER :: numpw!number of states (products of wanniers)
    INTEGER :: nums_psi!number of states (KS)
    REAL(kind=DP), DIMENSION(:,:,:), POINTER :: wwp!terms 
 END TYPE wp_psi
 
 TYPE wannier_u_prim
!this structure describes the  transformation
!from KS to Wannier states in the manifold C'
    INTEGER :: nums!total number of states
    INTEGER :: nums_occ!number of occupied states
    INTEGER :: nums_prim!number of states in manifold C'
    COMPLEX(kind=DP), DIMENSION(:,:), POINTER :: umat!INVERSE transformation matrix to wannier: Psi_c'=U_{c',j}w_j
 END TYPE wannier_u_prim

 TYPE v_pot_prim
!this structure describes the terms <\tilde{w}^P'_i|V|\tilde{w}^P_j'>
    INTEGER ::  numpw!number of states (products of wanniers)
    INTEGER ::  numpw_prim!number of states in manifold C' (products of wanniers)
    INTEGER, DIMENSION(:,:), POINTER :: ij!array of dimesion(2,numpw_prim) defining the product of w^C'*w^C
    REAL(kind=DP), DIMENSION(:,:), POINTER :: vmat!coulumbian matrix
    LOGICAL :: is_parallel!if true is a parallel part of the global cprim_prod matrix on polarization basis
    INTEGER :: numpw_para!numer of states(orthonormalized products) on this processor
    INTEGER :: first_para!first state (orthonormalized products, global order)on this processor
 END TYPE v_pot_prim

 TYPE wp_psi_cutoff_index
!this structure contains the indices for the description of terms \int Psi_i(r)\tilde{w}^P_i\tilde{w}^P_jPsi_i(r)dr
    INTEGER :: numpw!number of states (products of wanniers)
    INTEGER :: nums_psi!number of states (KS)
    INTEGER :: numpwpw!number of products of wannier products
    INTEGER, DIMENSION(:,:), POINTER :: index! of dimension (2,numpwpw) indices to wannier products
 END TYPE wp_psi_cutoff_index

 TYPE wp_psi_cutoff_data
!this structure contains the data for the description of terms \int Psi_i(r)\tilde{w}^P_i\tilde{w}^P_jPsi_i(r)dr
    INTEGER :: numpw!number of states (products of wanniers)
    INTEGER :: nums_psi!number of states (KS)
    INTEGER :: numpwpw!number of products of wannier products
    REAL(kind=DP), DIMENSION(:,:), POINTER :: wwp!  of dimension (numpwpw,nums_psi) 
 END TYPE wp_psi_cutoff_data

 TYPE  head_epsilon
!this structure contains the data for the descrpition of the head of the dielectric matrix
!calculated accurately with k_points sampling
!it also contains the data for the treatment of the wings
    INTEGER :: n!number of frequency steps
    REAL(kind=DP) :: omega!frequency range
    REAL(kind=DP), DIMENSION(:), POINTER :: freqs!frequency steps 2n+1
    REAL(kind=DP), DIMENSION(:,:), POINTER :: head!elements G=0,G=0 of the dielectric matrix
    INTEGER :: numpw!number of products of wanniers
    REAL(kind=DP), DIMENSION(:), POINTER :: gzero!G=0 elements of non orthogonal products of wanniers \tilde{w^P}
    REAL(kind=DP), DIMENSION(:,:,:), POINTER :: wing!contains the terms \Sum_G <G|\tilde{w^P_i}>\epsilon(G,G'=0; iw)
    REAL(kind=DP), DIMENSION(:,:,:), POINTER :: wing_c!contains the terms \Sum_G <G|\tilde{w^P_i}>\epsilon(G=0,G'; iw)
 END TYPE head_epsilon

 TYPE  cprim_prod
!this structure contains the terms \int Psi_c'(r) Psi_c(r) v(r,r') \tilde{w^P_i}dr dr'
!it can contain also the terms \int Psi_i(r) Psi_v,c(r) v(r,r') \tilde{w^P_i}dr dr'
!it can contain also the terms \int Psi_v(r) Psi_c(r) \tilde{w^P_i} dr
    INTEGER :: cprim!conduction band considered
    INTEGER :: nums!total number of states
    INTEGER :: nums_occ!number of occupied states
    INTEGER :: nums_cond!total number of conduction states
    INTEGER :: numpw!number of products of wanniers
    REAL(kind=DP), DIMENSION(:,:), POINTER :: cpmat!product terms
    INTEGER  :: lda!leading dimension of cpmat important for parallel execution
    LOGICAL :: is_parallel!if true is a parallel part of the global cprim_prod matrix on polarization basis
    INTEGER :: numpw_para!numer of states(orthonormalized products) on this processor
    INTEGER :: first_para!first state (orthonormalized products, global order)on this processor

 END TYPE cprim_prod


 TYPE upper_states
!this structure contains the data for the reduced upper states
    INTEGER :: nums!total number of REGULAR states 
    INTEGER :: nums_occ!number of occupied states
    INTEGER :: nums_reduced!number of reduced states 
    INTEGER :: nums_tot!number of TOTAL states
    REAL(kind=DP), DIMENSION(:), POINTER :: ene!KS energies of reduced states
 END TYPE upper_states

 TYPE vt_mat_lanczos
!this structure describes the terms   
!V^v_{v,l}=<Pc w_v(r)w^P_i(r)|z^v_l>
!where {z^v_l} is an orthonormal basis set which depends on v                                                               
    INTEGER :: ii!state v 
    INTEGER :: nums_occ!number of valence states 
    INTEGER :: numpw!dimension of polarization basis
    INTEGER :: numl!number orthonormal states {z^v_l}
    REAL(kind=DP), DIMENSION(:,:), POINTER :: vt_mat!matrix (numpw, numl)
 END TYPE vt_mat_lanczos

 TYPE tt_mat_lanczos
!this structure describes the terms 
!T^v_{i,j}=<z^v_i|t_j>
!where {t_j} is an orthonormal basis set spanning the whole manifold
!of the {z^v_l} for all the v
     INTEGER :: numt!dimension of the basis {t_j}
     INTEGER :: numl!number orthonormal states {z^v_l} 
     INTEGER :: ii!state v
     REAL(kind=DP), DIMENSION(:,:), POINTER :: tt_mat!matrix (numt,numl)
  END TYPE tt_mat_lanczos

  TYPE mat_lanczos_full
!this structures describes the terms
!M^{i,s}_{\mu\alpha}=<\psi_{i,s}(v\Phi_\mu)|\svev_\alpha}
     INTEGER :: ii!state KS
     INTEGER :: numpw!dimension of polarization basis        
     INTEGER :: nums!number of global s vectors
     COMPLEX(kind=DP), DIMENSION(:,:,:), POINTER :: f_mat(:,:,:)!(numpw,nums,2)
  END TYPE mat_lanczos_full
  TYPE lanczos_chain
!this structure described the lanczos chains and relative overlap
!starting from a basis set {t_j}
     INTEGER :: numt!dimension of the basis {t_j} 
     INTEGER :: num_steps!number of lanczos steps
     INTEGER :: ii!index of corresponding KS state not used for polarization
     REAL(kind=DP), DIMENSION(:,:,:), POINTER :: o_mat! (numt,num_steps,numt) overlaps <t_i|s^j_l>
                                                      !with s^j_l l-th lanczos vector staring from t_j
     REAL(kind=DP), DIMENSION(:,:), POINTER :: d!diagonal terms of H operator (num_steps,numt)
     REAL(kind=DP), DIMENSION(:,:), POINTER :: f!upper diagonal terms of H operator (num_steps,numt)
  END TYPE lanczos_chain

  TYPE partial_occ
!this structure described the date for treating partially occupied states when calculating P
     INTEGER :: nums_occ!total number of occupied states (also partially)
     INTEGER :: nums_occ_min!total number of fully occupied states
     INTEGER :: numpw!dimension of polarizability basis
     REAL(kind=DP), DIMENSION(:), POINTER :: f_occ!occupations of KS states
     REAL(kind=DP), DIMENSION(:,:,:), POINTER :: ppp_mat!overlaps psi_v(r)psi_v'(r)phi_mu(r)

  END TYPE partial_occ


  TYPE semicore
!this structure contains the terms \inr dr psi_i(r)\psi^sc_v(r)\Phi_mu(r)
!and the energies of semicore states
     INTEGER :: numpw!dimension of polarizability basis
     INTEGER :: n_semicore!number of semicore states
     INTEGER :: nums!number of KS states
     REAL(kind=DP), DIMENSION(:), POINTER :: en_sc!semicore energies in Ry
     REAL(kind=DP), DIMENSION(:,:,:), POINTER :: ppw_mat!overlaps dimension:numpw,n_semicore,nums

  END TYPE semicore


  TYPE contour_terms
!this structure contains the terms <psi_i|s_\alpha>
     INTEGER :: nums!number of KS states 
     INTEGER :: numt!dimension of global s basis
     REAL(kind=DP), DIMENSION(:,:),POINTER :: cmat!the overlaps
  END TYPE contour_terms

  TYPE full_prods
!this structure contains the terms \int dr conjg(\psi_i(r))\psi_j(r)(v\phi_\mu)(r)
!for full relativistic calculations 
     INTEGER :: nums!number of KS states of iterest
     INTEGER :: nbnd!total number of KS states
     INTEGER :: numpw!dimension of polarizability basis
     INTEGER :: numv!number of occupied valence states
     REAL(kind=DP), DIMENSION(:), POINTER :: ene_ks!KS energies (nbnd)
     COMPLEX(kind=DP), DIMENSION(:,:,:,:), POINTER :: gmat!product terms (numpw,2,nbnd,nums)
  END TYPE full_prods


!these routines deallocate the allocates structures
      INTERFACE free_memory

         MODULE PROCEDURE free_wannier_u,free_wannier_P,free_v_pot, free_q_mat,  free_memory_ortho_polaw, &
              &free_memory_wp_psi, free_memory_wannier_u_prim, free_memory_v_pot_prim, free_memory_wp_psi_cutoff_index,&
              &free_memory_wp_psi_cutoff_data, free_memory_head_epsilon, free_cprim_prod, free_memory_upper_states,&
              free_memory_vt_mat_lanczos, free_memory_tt_mat_lanczos, free_memory_lanczos_chain, free_memory_partial_occ,&
              free_memory_semicore,free_memory_contour_terms,free_memory_mat_lanczos_full,free_memory_full_prods
         
      END INTERFACE


      INTERFACE initialize_memory
         MODULE PROCEDURE initialize_memory_cprim_prod, initialize_memory_upper_states,initialize_memory_vt_mat_lanczos,&
              initialize_memory_tt_mat_lanczos,initialize_memory_lanczos_chain,initialize_memory_partial_occ,&
              initialize_memory_semicore,initialize_memory_contour_terms,initialize_memory_mat_lanczos_full,&
              initialize_memory_full_prods
      END INTERFACE


    



    CONTAINS

      subroutine free_wannier_u( r)
        implicit none
        type(wannier_u) :: r
        if(associated(r%ene))  deallocate(r%ene)
        nullify(r%ene)
        if(associated(r%umat)) deallocate(r%umat)
        nullify(r%umat)
        if(associated(r%ene_xc))  deallocate(r%ene_xc)
        nullify(r%ene_xc)
         if(associated(r%ene_lda_h))  deallocate(r%ene_lda_h)
        nullify(r%ene_lda_h)
        return
      end subroutine free_wannier_u


      subroutine free_wannier_P( w_P)
        implicit none
        type(wannier_P) :: w_P
        if(associated(w_P%ij)) deallocate(w_P%ij)
        nullify(w_P%ij)
        if(associated(w_P%o))  deallocate(w_P%o)
        nullify(w_P%o) 
        return
      end subroutine free_wannier_P

      subroutine free_v_pot(vp)
        implicit none
        type(v_pot) :: vp
        if(associated(vp%vmat)) deallocate(vp%vmat)
        nullify(vp%vmat)
        return
      end subroutine free_v_pot


      subroutine free_q_mat( qm)
        implicit none
        type(q_mat) :: qm
        integer :: iw
        if(associated(qm%wp)) then
           do iw=1,qm%numpw_para
              call free_wannier_P(qm%wp(iw))
           enddo
           deallocate(qm%wp)
           nullify(qm%wp)
        endif
        return
      end subroutine free_q_mat
      
      SUBROUTINE free_memory_ortho_polaw(op)
        !this subroutine deallocates the green descriptor
        implicit none
        TYPE(ortho_polaw) op
        if(associated(op%on_mat)) deallocate(op%on_mat)
        nullify(op%on_mat)
        return
      END SUBROUTINE free_memory_ortho_polaw

      SUBROUTINE free_memory_wp_psi( wp)
        implicit none
        TYPE(wp_psi) :: wp
        if(associated(wp%wwp))  deallocate(wp%wwp)
        nullify(wp%wwp)
        return
      END SUBROUTINE free_memory_wp_psi

      SUBROUTINE  free_memory_wannier_u_prim( ww)
        implicit none
        TYPE(wannier_u_prim) :: ww
        if(associated(ww%umat)) deallocate(ww%umat)
        nullify(ww%umat)
        return
      END SUBROUTINE free_memory_wannier_u_prim

      SUBROUTINE free_memory_v_pot_prim(vp)
        implicit none
        TYPE(v_pot_prim) :: vp

        if(associated(vp%ij)) deallocate(vp%ij)
        nullify(vp%ij)
        if(associated(vp%vmat)) deallocate(vp%vmat)
        nullify(vp%vmat)
      END SUBROUTINE free_memory_v_pot_prim
        
      SUBROUTINE free_memory_wp_psi_cutoff_index(wpi)
        implicit none
        TYPE(wp_psi_cutoff_index) :: wpi
        if(associated(wpi%index)) deallocate(wpi%index)
        nullify(wpi%index)
        return
      END SUBROUTINE free_memory_wp_psi_cutoff_index

      SUBROUTINE free_memory_wp_psi_cutoff_data(wp)
        implicit none
        TYPE(wp_psi_cutoff_data) :: wp
        if(associated(wp%wwp)) deallocate(wp%wwp)
        nullify(wp%wwp)
        return
      END SUBROUTINE free_memory_wp_psi_cutoff_data

      SUBROUTINE free_memory_head_epsilon(he)
        implicit none
        TYPE(head_epsilon) :: he
        if(associated(he%freqs)) deallocate(he%freqs)
        nullify(he%freqs)
        if(associated(he%head)) deallocate(he%head)
        nullify(he%head)
        if(associated(he%gzero)) deallocate(he%gzero)
        nullify(he%gzero)
        if(associated(he%wing)) deallocate(he%wing)
        nullify(he%wing)
        if (associated(he%wing_c)) deallocate(he%wing_c)
        nullify(he%wing_c)

        
      END SUBROUTINE free_memory_head_epsilon


      SUBROUTINE free_cprim_prod(cpp)
        implicit none
        TYPE(cprim_prod) :: cpp

        if(associated(cpp%cpmat)) deallocate( cpp%cpmat)
        nullify(cpp%cpmat)

        return
      END SUBROUTINE free_cprim_prod

      SUBROUTINE free_memory_full_prods(fp)
        implicit none
        TYPE(full_prods) ::fp
        if(associated(fp%ene_ks)) deallocate(fp%ene_ks)
        nullify(fp%ene_ks)
        if(associated(fp%gmat)) deallocate(fp%gmat)
        nullify(fp%gmat)

        return
      END SUBROUTINE free_memory_full_prods

      SUBROUTINE initialize_memory_full_prods(fp)
        implicit none
        TYPE(full_prods) ::fp
        nullify(fp%ene_ks)
        nullify(fp%gmat)

        return
      END SUBROUTINE initialize_memory_full_prods



      SUBROUTINE initialize_memory_cprim_prod(cpp)
         implicit none
        TYPE(cprim_prod) :: cpp

        nullify(cpp%cpmat)

        return
      END SUBROUTINE initialize_memory_cprim_prod

      SUBROUTINE free_memory_upper_states(us)
        implicit none
        TYPE(upper_states) :: us
        if(associated(us%ene)) deallocate(us%ene)
        nullify(us%ene)
      END SUBROUTINE free_memory_upper_states



      SUBROUTINE free_memory_partial_occ(po)
        implicit none
        TYPE(partial_occ) :: po
        if(associated(po%f_occ)) deallocate(po%f_occ)
        nullify(po%f_occ)
        if(associated(po%ppp_mat)) deallocate(po%ppp_mat)
        nullify(po%ppp_mat)
      END SUBROUTINE free_memory_partial_occ

      SUBROUTINE free_memory_semicore(sc)
        implicit none
        TYPE(semicore) :: sc
        if(associated(sc%en_sc)) deallocate(sc%en_sc)
        nullify(sc%en_sc)
        if(associated(sc%ppw_mat)) deallocate(sc%ppw_mat)
        nullify(sc%ppw_mat)
      END SUBROUTINE free_memory_semicore


      SUBROUTINE initialize_memory_semicore(sc)
        implicit none
        TYPE(semicore) :: sc
        nullify(sc%en_sc)
        nullify(sc%ppw_mat)
      END SUBROUTINE initialize_memory_semicore


      SUBROUTINE free_memory_contour_terms(ct)
        implicit none
        TYPE(contour_terms) :: ct
        if(associated(ct%cmat)) deallocate(ct%cmat)
        nullify(ct%cmat)
      END SUBROUTINE free_memory_contour_terms

      SUBROUTINE initialize_memory_contour_terms(ct)
        implicit none
        TYPE(contour_terms) :: ct
        nullify(ct%cmat)
      END SUBROUTINE initialize_memory_contour_terms

      SUBROUTINE initialize_memory_partial_occ(po)
        implicit none
        TYPE(partial_occ) :: po
        nullify(po%f_occ)
        nullify(po%ppp_mat)
      END SUBROUTINE initialize_memory_partial_occ



       SUBROUTINE initialize_memory_upper_states(us)
        implicit none
        TYPE(upper_states) :: us
        nullify(us%ene)
      END SUBROUTINE initialize_memory_upper_states
      
      SUBROUTINE free_memory_vt_mat_lanczos( vtl)
        implicit none
        TYPE(vt_mat_lanczos) :: vtl
        if(associated(vtl%vt_mat)) deallocate(vtl%vt_mat)
        nullify(vtl%vt_mat)
      END SUBROUTINE free_memory_vt_mat_lanczos

      SUBROUTINE initialize_memory_vt_mat_lanczos( vtl)
         implicit none
         TYPE(vt_mat_lanczos) :: vtl
         nullify(vtl%vt_mat)
       END SUBROUTINE initialize_memory_vt_mat_lanczos


       SUBROUTINE free_memory_mat_lanczos_full( full)
         implicit none
        TYPE(mat_lanczos_full) :: full
        if(associated(full%f_mat)) deallocate(full%f_mat)
        nullify(full%f_mat)
      END SUBROUTINE free_memory_mat_lanczos_full
      
      SUBROUTINE initialize_memory_mat_lanczos_full( full)
        implicit none
        TYPE(mat_lanczos_full) :: full
        nullify(full%f_mat)
      END SUBROUTINE initialize_memory_mat_lanczos_full


       SUBROUTINE free_memory_tt_mat_lanczos( ttl)
        implicit none
        TYPE(tt_mat_lanczos) :: ttl
        if(associated(ttl%tt_mat)) deallocate(ttl%tt_mat)
        nullify(ttl%tt_mat)
      END SUBROUTINE free_memory_tt_mat_lanczos

      SUBROUTINE initialize_memory_tt_mat_lanczos( ttl)
         implicit none
         TYPE(tt_mat_lanczos) :: ttl
         nullify(ttl%tt_mat)
       END SUBROUTINE initialize_memory_tt_mat_lanczos

       SUBROUTINE initialize_memory_lanczos_chain( lc)
          implicit none
         TYPE(lanczos_chain) :: lc
         nullify(lc%o_mat)
         nullify(lc%d)
         nullify(lc%f)
       END SUBROUTINE initialize_memory_lanczos_chain

       SUBROUTINE free_memory_lanczos_chain( lc)
          implicit none
         TYPE(lanczos_chain) :: lc
         if(associated(lc%o_mat)) deallocate(lc%o_mat)
         if(associated(lc%d)) deallocate(lc%d)
         if(associated(lc%f)) deallocate(lc%f)
         nullify(lc%o_mat)
         nullify(lc%d)
         nullify(lc%f)
       END SUBROUTINE free_memory_lanczos_chain

       
 END MODULE basic_structures



